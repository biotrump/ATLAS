#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_level1.h"
#include "atlas_tlvl3.h"
void Mjoin(PATL,DoWork_amm_tMN)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tamm_tMN_t *pd = pp->PD;
   cm2am_t a2blk = pd->a2blk, b2blk = pd->b2blk;
   ammkern_t amm = pd->amm_b0;
   TYPE *wA = pd->w + vrank*pd->szW;
   TYPE *wB;
   TYPE *wC;
   const int nkblks = pd->nkblks, kb=pd->kb, KB0=pd->KB0, kb0=pd->kb0,
             M=pd->M, N=pd->N;
   const size_t incA = (pd->TA) ? 1:pd->lda;
   const size_t incB = (pd->TB) ? pd->ldb:1;
   int kctr;

   wA = ATL_AlignPtr(wA);
   wB = wA + pd->szA;
   wC = wB + pd->szB;
   while ((kctr = ATL_DecGlobalAtomicCount(pd->KbCtr, vrank)))
   {
      const int kblk = nkblks - kctr;
      size_t k;
/*
 *    Normal case, doing full K blocks
 */
      if (kblk)
      {
         size_t k = kb0 + kb*(kblk-1);
         const TYPE *a = pd->A + k*incA, *b = pd->B + k*incB;
         a2blk(kb, M, ATL_rone, a, pd->lda, wA);
         b2blk(kb, N, ATL_rone, b, pd->ldb, wB);
         amm(pd->nmu, pd->nnu, kb, wA, wB, wC, wA, wB, wC);
      }
/*
 *    Doing possibly partial kb0 at beginning
 */
      else
      {
         amm = (amm == pd->amm_b1) ? pd->ammK_b1 : pd->ammK_b0;
         a2blk(kb0, M, ATL_rone, pd->A, pd->lda, wA);
         b2blk(kb0, N, ATL_rone, pd->B, pd->ldb, wB);
         amm(pd->nmu, pd->nnu, KB0, wA, wB, wC, wA, wB, wC);
      }
      amm = pd->amm_b1;
   }
/*
 * If I did no work, zero my wC so combine isn't messed up!
 */
   if (amm == pd->amm_b0)
      Mjoin(PATL,zero)(pd->nC, wC, 1);
}

void Mjoin(PATL,DoComb_amm_tMN)
   (void *vpp, int rank, int vrank, int hisvrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tamm_tMN_t *pd = pp->PD;
   TYPE *mC, *hC;

   mC = pd->w + vrank*pd->szW;
   mC = ATL_AlignPtr(mC);
   mC += pd->szA + pd->szB;

   hC = pd->w + hisvrank*pd->szW;
   hC = ATL_AlignPtr(hC);
   hC += pd->szA + pd->szB;
   Mjoin(PATL,axpy)(pd->nC, ATL_rone, hC, 1, mC, 1);
}

/*
 * This routine handles the case where M <= maxMB && N <= maxNB so C is
 * only one block.  It gets its own case because it requires minimal
 * workspace.  K will need to fairly long, and due to streaming A & B
 * through the cache with only intra-kernel reuse, its parallel efficiency
 * tends to be limited by bus speed.
 */
int Mjoin(PATL,tammm_tMN)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CINT M,
   ATL_CINT N,
   ATL_CINT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CINT lda,
   const TYPE *B,
   ATL_CINT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CINT ldc
)
{
   ATL_SZT nmblks;
   amminfo_t mminfo;
   unsigned int i, mb, nb, kb, kb0, KB0, nmu, nnu, mu, nu, ku, P, mr;
   ATL_INT nkblks;
   ATL_tamm_tMN_t pd;  /* problem definition structure */
   void *vp;
   TYPE *wC;

/*
 * Special case for tiny N&K, and large M
 */
   if (N > ATL_AMM_MAXNB || M > ATL_AMM_MAXMB ||
       K < Mmin(8,ATL_NTHREADS)*ATL_AMM_MAXKB)
      return(1);
/*
 * Not uncommon to recur and hit this degenerate case that's a dot product
 */
   if (M == 1 && N == 1)
   {
      TYPE dot;
      dot = Mjoin(PATL,dot)(K, A, (TA == AtlasNoTrans)?lda:1,
                            B, (TB == AtlasNoTrans) ? 1:ldb);
      if (SCALAR_IS_ZERO(beta))
         *C = alpha * dot;
      else
         *C = alpha * dot + beta * (*C);
      return(0);
   }
   nkblks = K / ATL_AMM_66KB;
   mb = Mjoin(PATL,tGetAmmmInfo)(&mminfo, Mmin(nkblks, ATL_NTHREADS), TA,
                                 TB, M, N, K, alpha, beta);

   if (!SCALAR_IS_ONE(alpha)) { ATL_assert(mb == 2); }
   pd.a2blk = mminfo.a2blk;
   pd.b2blk = mminfo.b2blk;
   pd.amm_b0 = mminfo.amm_b0;
   pd.amm_b1 = mminfo.amm_b1;
   pd.TA = (TA == AtlasTrans);
   pd.TB = (TB == AtlasTrans);
   pd.M = M;
   pd.N = N;
   pd.K = K;
   pd.A = A;
   pd.B = B;
   pd.C = C;
   pd.lda = lda;
   pd.ldb = ldb;
   pd.ldc = ldc;
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   pd.nmu = nmu = (M+mu-1)/mu;
   pd.nnu = nnu = (N+nu-1)/nu;
   kb = pd.kb = mminfo.kb;
   nkblks = K / kb;
   KB0 = kb0 = K - nkblks*kb;
   if (!kb0)
   {
      kb0 = KB0 = kb;
      pd.ammK_b0 = mminfo.amm_b0;
      pd.ammK_b1 = mminfo.amm_b1;
   }
   else
   {
      nkblks++;
      #if ATL_AMM_MAXKMAJ > 1
         if (ATL_AMMFLG_KMAJOR(mminfo.flag))
         {
            KB0 = ((kb0+ku-1)/ku)*ku;
            if (ATL_AMMFLG_KRUNTIME(mminfo.flag))
            {
               pd.ammK_b0 = mminfo.amm_b0;
               pd.ammK_b1 = mminfo.amm_b1;
            }
            else
            {
               pd.ammK_b0 = mminfo.amm_k1_b0;
               pd.ammK_b1 = mminfo.amm_k1_b1;
            }
         }
         else
      #endif
      {
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag) && kb0 == (kb0/ku)*ku &&
             kb0 > mminfo.kbmin)
         {
            pd.ammK_b0 = mminfo.amm_b0;
            pd.ammK_b1 = mminfo.amm_b1;
         }
         else
         {
            pd.ammK_b0 = mminfo.amm_k1_b0;
            pd.ammK_b1 = mminfo.amm_k1_b1;
         }
      }
   }
   pd.nkblks = nkblks;
   pd.KB0 = KB0;
   pd.kb0 = kb0;
/*
 * Maximum scale is limited by NTHREADS or max number of K-blocks
 */
   P = (ATL_NTHREADS <= nkblks) ? ATL_NTHREADS : nkblks;
   mb = nmu * mu;
   nb = nnu * nu;
   pd.nC = mb*nb;
   pd.szA = mb*kb;
   pd.szB = nb*kb;
   pd.szW = pd.szA + pd.szB + pd.nC + ATL_Cachelen;

   vp = malloc(ATL_MulBySize(pd.szW*P + mu*nu*ku));
   if (!vp)
      return(2);
   pd.w = vp;
   pd.KbCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nkblks, P), nkblks, 0);

/*   #define DEBUG1 */
   #ifdef DEBUG1
   {
      ATL_tpool_t *pp=ATL_TP_PTR;
      if (!pp)
         pp = ATL_NewThreadPool(1, 0, NULL);
      ATL_assert(pp);
      pp->PD = &pd;
      Mjoin(PATL,DoWork_amm_tMN)(pp, 0, 0);
      if (pp != ATL_TP_PTR)
         ATL_FreeThreadPool(pp);
   }
   #else
      ATL_goParallel(P, Mjoin(PATL,DoWork_amm_tMN), Mjoin(PATL,DoComb_amm_tMN),
                     &pd, NULL);
   #endif
   ATL_FreeGlobalAtomicCount(pd.KbCtr);
/*
 * Copy answer back out while scaling by alpha and beta
 */
   wC = ATL_AlignPtr(pd.w);
   wC += pd.szA + pd.szB;
   mminfo.Cblk2cm(M, N, alpha, wC, beta, C, ldc);

   free(vp);
   return(0);
}
