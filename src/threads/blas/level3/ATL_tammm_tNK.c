#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"

void Mjoin(PATL,DoWork_tamm_tNK)(ATL_LAUNCHSTRUCT_t *lp, void *vp)
{
   ATL_thread_t *tp = vp;
   ATL_tamm_tNK_t *pd = lp->opstruct;  /* problem definition structure */
   const unsigned int rank = tp->rank, mb=pd->mb, kb=pd->KB0, K=pd->K,
      N=pd->N, nmblks = pd->nmblks, nmbm1=nmblks-1, nmu=pd->nmu, nnu=pd->nnu,
      lda=pd->lda, ldc=pd->ldc;
   TYPE *pB, *pA, *pC, *C = pd->C;
   const TYPE *A = pd->A;
   int BCOPIED=0;
   int imtr;
   const TYPE beta = *pd->beta;
   const size_t incA = (pd->TA) ? lda : 1;
   cm2am_t a2blk = pd->a2blk;
   ablk2cmat_t blk2c = pd->blk2c;
   ammkern_t amm = pd->amm_b0;

   pB = pd->w;
/*
 * First guy here starts to copy B
 */
   if (ATL_DecAtomicCount(pd->BassgCtr))
   {
      pd->b2blk(K, N, *pd->alpha, pd->B, pd->ldb, pB);
/*
 *    Let waiting threads know B is ready for use
 */
      BCOPIED = ATL_DecAtomicCount(pd->BdoneCtr);
      ATL_assert(BCOPIED);
   }
   pA = pB + pd->bsz + rank*(pd->wsz);
   pA = ATL_AlignPtr(pA);
   pC = pA + mb*kb;
   pC = ATL_AlignPtr(pC);
/*
 * For first block I work on, I must await B to be copied
 */
   if (!BCOPIED)
   {
      int iblk;
      size_t ii;

      imtr = ATL_DecGlobalAtomicCount(pd->MbCtr, rank);
      if (imtr)
      {
         iblk = nmblks - imtr;
         ii = mb*iblk;
         if (iblk != nmbm1)
            a2blk(K, mb, ATL_rone, A+incA*ii, lda, pA);
         else
            a2blk(K, pd->mr, ATL_rone, A+incA*ii, lda, pA);
         while (ATL_GetAtomicCount(pd->BdoneCtr))  /* await B cpy finish */
            ATL_thread_yield();
         if (iblk != nmbm1)
         {
            amm(nmu, nnu, kb, pA, pB, pC, pA, pB, pC);
            blk2c(mb, N, ATL_rone, pC, beta, C+ii, ldc);
         }
         else
         {
            amm(pd->nmuL, nnu, kb, pA, pB, pC, pA, pB, pC);
            blk2c(pd->mr, N, ATL_rone, pC, beta, C+ii, ldc);
         }
      }
   }
/*
 * Now, B is ready, so just go to town on remaining blocks
 */
   while ((imtr = ATL_DecGlobalAtomicCount(pd->MbCtr, rank)))
   {
      const int iblk = nmblks - imtr;
      const size_t ii = mb*iblk;

      if (iblk != nmbm1)
      {
         a2blk(K, mb, ATL_rone, A+incA*ii, lda, pA);
         amm(nmu, nnu, kb, pA, pB, pC, pA, pB, pC);
         blk2c(mb, N, ATL_rone, pC, beta, C+ii, ldc);
      }
      else
      {
         a2blk(K, pd->mr, ATL_rone, A+incA*ii, lda, pA);
         amm(pd->nmuL, nnu, kb, pA, pB, pC, pA, pB, pC);
         blk2c(pd->mr, N, ATL_rone, pC, beta, C+ii, ldc);
      }
   }
}
/*
 * This routine handles the case where N <= maxNB && K <= maxKB, so B is
 * only one block.  It is particularly important for the panel factorizations
 * of both LU and QR.
 */
int Mjoin(PATL,tammm_tNK)
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
   unsigned int i, mb, nb, kb, mu, nu, ku, P, mr;
   ATL_tamm_tNK_t pd;  /* problem definition structure */
   void *vp;

/*
 * Special case for tiny N&K, and large M
 */
   if (N >= ATL_AMM_MAXNB || K >= ATL_AMM_MAXKB || M < ATL_AMM_MAXMB ||
       M < Mmin(8,ATL_NTHREADS)*ATL_AMM_MAXMB)
      return(1);
   Mjoin(PATL,GetRankKInfo)(&mminfo, TA, TB, M, N, K, alpha, beta);
   pd.a2blk = mminfo.a2blk;
   pd.b2blk = mminfo.b2blk;
   pd.blk2c = mminfo.Cblk2cm;
   pd.amm_b0 = mminfo.amm_b0;
   pd.TA = (TA == AtlasTrans);
   pd.TB = (TB == AtlasTrans);
   pd.N = N;
   pd.K = K;
   pd.A = A;
   pd.B = B;
   pd.C = C;
   pd.lda = lda;
   pd.ldb = ldb;
   pd.ldc = ldc;
   pd.alpha = &alpha;
   pd.beta  = &beta;
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   pd.mb = mb = mminfo.mb;
   pd.nmu = mb / mu;
   pd.nnu = (N+nu-1)/nu;
   nb = pd.nnu * nu;
   kb = mminfo.kb;
   nmblks = M / mb;
   mr = M - nmblks*mb;
   if (!mr)
   {
      pd.mbL = mr = mb;
      pd.nmuL = pd.nmu;
   }
   else
   {
      nmblks++;
      pd.nmuL = (mr+mu-1)/mu;
      pd.mbL = pd.nmuL * mu;
   }
   pd.mr = mr;
   pd.nmblks = nmblks;
   pd.KB0 = K;
   #if ATL_MAXKMAJ_RKK > 1
      if (ATL_AMMFLG_KMAJOR(mminfo.flag))
         pd.KB0 = ((K+ku-1)/ku)*ku;
   #endif
/*
 * Maximum scale is limited by NTHREADS or max number of M-blocks
 */
   P = (ATL_NTHREADS <= nmblks) ? ATL_NTHREADS : nmblks;
/*
 * We have a common B wrk of size KB0*nb, then
 * for each node, we need workspace: sz(A,C) = mb*K, K*nb, mb*N, laid out
 * in memory as A,C, then we add safety margin mu*nu*ku so advance loads don't
 * seg fault, and we add space for aligning the ptrs
 */
   pd.bsz = pd.KB0*nb;
   pd.wsz = mb*(pd.nnu*nu + pd.bsz) + 2*ATL_DivBySize(ATL_Cachelen);

   vp = malloc(ATL_MulBySize(pd.wsz*P + pd.bsz+mu*nu*ku) + ATL_Cachelen);
   if (!vp)
      return(2);
   pd.w = ATL_AlignPtr(vp);
   pd.MbCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nmblks, P), nmblks, 0);
   pd.BassgCtr = ATL_SetAtomicCount(1);
   pd.BdoneCtr = ATL_SetAtomicCount(1);
   #ifdef DEBUG1
   {
      ATL_LAUNCHSTRUCT_t ls;
      ATL_thread_t ts;
      ts.rank = 0;
      ts.P = 1;
      ls.opstruct = &pd;
      Mjoin(PATL,DoWork_tamm_tNK)(&ls, &ts);
   }
   #else
      ATL_goparallel(P, Mjoin(PATL,DoWork_tamm_tNK), &pd, NULL);
   #endif

   ATL_FreeAtomicCount(pd.BdoneCtr);
   ATL_FreeAtomicCount(pd.BassgCtr);
   ATL_FreeGlobalAtomicCount(pd.MbCtr);

   free(vp);
   return(0);
}
