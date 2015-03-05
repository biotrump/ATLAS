#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_sum.h))

typedef struct ammrec ammrec_t;
struct ammrec
{
   cm2am_t a2blk, b2blk;
   ablk2cmat_t blk2c;
   ammkern_t amm_b0, amm_b1, amm_k1_b0, amm_k1_b1;
   size_t lda, ldb, ldc, incAm, incAk, incBk, incBn, incCn, incCm;
   TYPE alpA, alpB, alpC, beta;
   int mbkb, kbnb, mbnb, nmu, nnu, nmuF, nnuF;
   int mb, nb, kb, KB0, mr, nr, kr;
};


#define v_kp 1
#define v_np 2
#define v_mp 4
#define v_cpA 8
#define v_cpB 16
#define v_cpC 32
#define v_lwC 64  /* last guy to update C must write it out */
#define v_nmA 128 /* don't move A ptr */
#define v_nmB 256 /* don't move B ptr */
#define v_nmC 512 /* don't move C ptr */

#define b_kp 0
#define b_np 1
#define b_mp 2
#define b_cpA 3
#define b_cpB 4
#define b_cpC 5
#define b_lwC 6
#define b_nmA 7
#define b_nmB 8
#define b_nmC 9


static void ammmRECf /* no partial blocks */
(
   const ammrec_t *pd,
   int nmblks,
   int nnblks,
   int nkblks,
   int flag, /* bits: 0:kp, 1:np, 2:mp, 3:cpA, 4:cpyB, 5:cpyC */
   const TYPE *A,
   const TYPE *B,
   TYPE *C,
   TYPE *a,
   TYPE *b,
   TYPE *c
)
{
/*
 * Stop recursion and do multiply when only 1 block is left
 */
   if (nmblks == 1 && nnblks == 1 && nkblks == 1)
   {
      TYPE *an = (flag & v_nmA) ? a : a+pd->mbkb;
      TYPE *bn = (flag & v_nmB) ? b : b+pd->kbnb;
      TYPE *cn = (flag & v_nmC) ? c : c+pd->mbnb;
      const ammkern_t amm = (flag & v_cpC) ? pd->amm_b0 : pd->amm_b1;
      if (flag & v_cpA)
         pd->a2blk(pd->kb, pd->mb, pd->alpA, A, pd->lda, a);
      if (flag & v_cpB)
         pd->b2blk(pd->kb, pd->nb, pd->alpB, B, pd->ldb, b);
      amm(pd->nmu, pd->nnu, pd->kb, a, b, c, an, bn, cn);
      if (flag & v_lwC)
         pd->blk2c(pd->mb, pd->nb, pd->alpC, c, pd->beta, C, pd->ldc);
   }
   else if (nnblks >= nkblks && nnblks >= nmblks)   /* recursively divide N */
   {
      const int nR=(nnblks>>1), nL = nnblks-nR;
      ammmRECf(pd, nmblks, nL, nkblks, (flag|v_nmA) & ~(v_nmB+v_nmC),
               A, B, C, a, b, c);
      ammmRECf(pd, nmblks, nR, nkblks, flag & ~v_cpA,
              A, B+nL*pd->incBn, C+pd->incCn*nL,
              a, b+nL*nkblks*pd->kbnb, c+nL*nmblks*pd->mbnb);
   }
   else if (nmblks >= nkblks)                       /* recursively divide M */
   {
      const int nR=(nmblks>>1), nL = nmblks-nR;
      ammmRECf(pd, nL, nnblks, nkblks, (flag|v_nmB) & ~(v_nmA+v_nmC),
               A, B, C, a, b, c);
      ammmRECf(pd, nR, nnblks, nkblks, flag & ~v_cpB,
              A+nL*pd->incAm, B, C+pd->mb*nL,
              a+nL*nkblks*pd->mbkb, b, c+nL*nnblks*pd->mbnb);
   }
   else                                             /* recursively divide K */
   {
      const int nR=(nkblks>>1), nL = nkblks-nR;
      ammmRECf(pd, nmblks, nnblks, nL, (flag|v_nmC) & ~(v_nmA+v_nmB+v_lwC),
               A, B, C, a, b, c);
      ammmRECf(pd, nmblks, nnblks, nR, flag & ~v_cpC,
               A+nL*pd->incAk, B+nL*pd->incBk, C,
               a+nL*nmblks*pd->mbkb, b+nL*nnblks*pd->kbnb, c);
   }
}

static void ammmREC  /* partial blocks are possible */
(
   const ammrec_t *pd,
   int nmblks,
   int nnblks,
   int nkblks,
   int flag, /* bits: 0:kp, 1:np, 2:mp, 3:cpA, 4:cpyB, 5:cpyC */
   const TYPE *A,
   const TYPE *B,
   TYPE *C,
   TYPE *a,
   TYPE *b,
   TYPE *c
)
{
/*
 * If I only have one possibly partial block left.  Equivalent to:
 * if ( ((nmblks == 1 && !(flag&v_mp)) || nmblks == 0) &&
 *      ((nnblks == 1 && !(flag&v_np)) || nnblks == 0) &&
 *      ((nkblks == 1 && !(flag&v_kp)) || nkblks == 0) )
 */
   if (nmblks < 2 && nnblks < 2 && nkblks < 2 &&
       nmblks != ((flag>>b_mp)&1) &&
       nnblks != ((flag>>b_np)&1) &&
       nkblks != ((flag>>b_kp)&1))
   {
      TYPE *an = (flag & v_nmA) ? a : a+pd->mbkb;
      TYPE *bn = (flag & v_nmB) ? b : b+pd->kbnb;
      TYPE *cn = (flag & v_nmC) ? c : c+pd->mbnb;
/*
 *    If we aren't doing K-cleanup, can always use fastest kernels
 */
      if (!(flag & v_kp))
      {
         const ammkern_t amm = (flag & v_cpC) ? pd->amm_b0 : pd->amm_b1;
         if (!(flag & (v_mp+v_np)))   /* one full block to multiply */
         {
            if (flag & v_cpA)
               pd->a2blk(pd->kb, pd->mb, pd->alpA, A, pd->lda, a);
            if (flag & v_cpB)
               pd->b2blk(pd->kb, pd->nb, pd->alpB, B, pd->ldb, b);
            amm(pd->nmu, pd->nnu, pd->kb, a, b, c, an, bn, cn);
            if (flag & v_lwC)
               pd->blk2c(pd->mb, pd->nb, pd->alpC, c, pd->beta, C, pd->ldc);
         }
         else
         {
            int nmu=pd->nmu, nnu=pd->nnu, m=pd->mb, n=pd->nb, kb=pd->kb;
            if (flag & v_mp)
            {
               nmu = pd->nmuF;
               m = pd->mr;
            }
            if (flag & v_np)
            {
               nnu = pd->nnuF;
               n = pd->nr;
            }
            if (flag & v_cpA)
               pd->a2blk(kb, m, pd->alpA, A, pd->lda, a);
            if (flag & v_cpB)
               pd->b2blk(kb, n, pd->alpB, B, pd->ldb, b);
            amm(nmu, nnu, kb, a, b, c, an, bn, cn);
            if (flag & v_lwC)
               pd->blk2c(m, n, pd->alpC, c, pd->beta, C, pd->ldc);
         }
      }
      else /* if (flag & (v_mp|v_np|v_kp))   partial block to multiply */
      {
         const ammkern_t amm = (flag & v_cpC) ? pd->amm_k1_b0 : pd->amm_k1_b1;
         const int nmu = (nmblks) ? pd->nmu : pd->nmuF;
         const int m = (nmblks) ? pd->mb : pd->mr;
         const int nnu = (nnblks) ? pd->nnu : pd->nnuF;
         const int n = (nnblks) ? pd->nb : pd->nr;
         const int k = (nkblks) ? pd->kb : pd->kr;
         const int KBF = (nkblks) ? k : pd->KB0;

         if (flag & v_cpA)
            pd->a2blk(k, m, pd->alpA, A, pd->lda, a);
         if (flag & v_cpB)
            pd->b2blk(k, n, pd->alpB, B, pd->ldb, b);
         amm(nmu, nnu, KBF, a, b, c, an, bn, cn);
         if (flag & v_lwC)
            pd->blk2c(m, n, pd->alpC, c, pd->beta, C, pd->ldc);
      }
   }
   else if (nnblks >= nkblks && nnblks >= nmblks &&           /* recursively */
            (nnblks > 1 || (flag&v_np)))                      /* divide N    */
   {
      const int nR=(nnblks>>1), nL = nnblks-nR;
      const size_t mblks=(flag&v_mp) ? nmblks+1 : nmblks;
      const size_t kblks=(flag&v_kp) ? nkblks+1 : nkblks;
      const int flg = (flag|v_nmA) & ~(v_np+v_nmB+v_nmC);
      if (flag & (v_mp+v_kp))
         ammmREC(pd, nmblks, nL, nkblks, flg, A, B, C, a, b, c);
      else
         ammmRECf(pd, nmblks, nL, nkblks, flg, A, B, C, a, b, c);
      ammmREC(pd, nmblks, nR, nkblks, flag & ~v_cpA,
              A, B+nL*pd->incBn, C+pd->incCn*nL,
              a, b+nL*kblks*pd->kbnb, c+nL*mblks*pd->mbnb);
   }
   else if (nmblks >= nkblks && (nmblks > 1 || (flag&v_mp)))  /* recursively */
   {                                                          /* divide M    */
      const int nR=(nmblks>>1), nL = nmblks-nR;
      const size_t nblks=(flag&v_np) ? nnblks+1 : nnblks;
      const size_t kblks=(flag&v_kp) ? nkblks+1 : nkblks;
      const int flg = (flag|v_nmB) & ~(v_mp+v_nmA+v_nmC);
      if (flag&(v_np+v_kp))
         ammmREC(pd, nL, nnblks, nkblks, flg, A, B, C, a, b, c);
      else
         ammmRECf(pd, nL, nnblks, nkblks, flg, A, B, C, a, b, c);
      ammmREC(pd, nR, nnblks, nkblks, flag & ~v_cpB,
              A+nL*pd->incAm, B, C+pd->mb*nL,
              a+nL*kblks*pd->mbkb, b, c+nL*nblks*pd->mbnb);
   }
   else                                                       /* recursively */
   {                                                          /* divide K    */
      const int nR=(nkblks>>1), nL = nkblks-nR;
      const size_t mblks=(flag&v_mp) ? nmblks+1 : nmblks;
      const size_t nblks=(flag&v_np) ? nnblks+1 : nnblks;
      const int flg = (flag|v_nmC) & ~(v_kp+v_lwC+v_nmA+v_nmB);
      if (flag & (v_mp+v_np))
         ammmREC(pd, nmblks, nnblks, nL, flg, A, B, C, a, b, c);
      else
         ammmRECf(pd, nmblks, nnblks, nL, flg, A, B, C, a, b, c);
      ammmREC(pd, nmblks, nnblks, nR, flag & ~v_cpC,
              A+nL*pd->incAk, B+nL*pd->incBk, C,
              a+nL*mblks*pd->mbkb, b+nL*nblks*pd->kbnb, c);
   }
}

int Mjoin(PATL,ammmREC)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   const TYPE *B,
   ATL_CSZT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   void *vp;
   TYPE *a, *b, *c;
   amminfo_t mminfo;
   ammrec_t pd;
   #if ATL_AMM_MAXKMAJ > 1
      size_t kb0U, KK;
   #else
      #define kb0U kb0
      #define KK K
   #endif
   size_t nmblks, nnblks, nkblks, szA, szB, szC, i, j, k;
   int mb, nb, kb, mu, nu, ku, mr, nr, kr, KB0, flag, appAl;

   pd.alpA = pd.alpB = pd.alpC = ATL_rone;
   appAl = Mjoin(PATL,GetAmmmInfo)(&mminfo, TA, TB, M, N, K, alpha, beta);
   pd.mb = mb = mminfo.mb;
   pd.nb = nb = mminfo.nb;
   pd.kb = kb = mminfo.kb;
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   nmblks = M / mb;
   nnblks = N / nb;
   nkblks = K / kb;
   pd.nmu = mb / mu;
   pd.nnu = mb / nu;
   pd.mr = mr = M - nmblks*mb;
   if (mr)
      pd.nmuF = (mr+mu-1)/mu;
   else
      pd.nmuF = pd.nmu;
   pd.nr = nr = N - nnblks*nb;
   if (nr)
      pd.nnuF = (nr+nu-1)/nu;
   else
      pd.nnuF = pd.nnu;
   pd.kr = kr = K - nkblks*kb;
   if (!appAl)
      pd.alpA = alpha;
   else if (appAl == 1)
      pd.alpB = alpha;
   else
      pd.alpC = alpha;
   pd.a2blk = mminfo.a2blk;
   pd.b2blk = mminfo.b2blk;
   pd.blk2c = mminfo.Cblk2cm;
   pd.amm_b0 = mminfo.amm_b0;
   pd.amm_b1 = mminfo.amm_b1;
   if (!kr)
   {
      pd.amm_k1_b0 = mminfo.amm_b0;
      pd.amm_k1_b1 = mminfo.amm_b1;
      KB0 = kb;
   }
/*
 * If last K-block is partial, compute K of gemm (KB0) and K of copy (kr)
 */
   else
   {
/*
 *    K-major require GEMM's K (KB0) to be a mult of ku
 */
      #if ATL_AMM_MAXKMAJ > 1
      if (ATL_AMMFLG_KMAJOR(mminfo.flag))
      {
         KB0 = ((kr+ku-1)/ku)*ku;
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag))
         {
            pd.amm_k1_b0 = mminfo.amm_b0;
            pd.amm_k1_b1 = mminfo.amm_b1;
         }
         else
         {
            pd.amm_k1_b0 = mminfo.amm_k1_b0;
            pd.amm_k1_b1 = mminfo.amm_k1_b1;
         }

      }
      else
      {
      #endif
         KB0 = kr;
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag) && kr == (kr/ku)*ku &&
             kr > mminfo.kbmin)
         {
            pd.amm_k1_b0 = mminfo.amm_b0;
            pd.amm_k1_b1 = mminfo.amm_b1;
         }
         else
         {
            pd.amm_k1_b0 = mminfo.amm_k1_b0;
            pd.amm_k1_b1 = mminfo.amm_k1_b1;
         }
      #if ATL_AMM_MAXKMAJ > 1
      }
      #endif
   }
   pd.mb = mminfo.mb;
   pd.nb = mminfo.nb;
   pd.kb = mminfo.kb;
   pd.mbkb = mminfo.mb * mminfo.kb;
   pd.kbnb = mminfo.kb * mminfo.nb;
   pd.mbnb = mminfo.mb * mminfo.nb;
   pd.lda = lda;
   pd.ldb = ldb;
   pd.ldc = ldc;
   pd.nmu = mminfo.mb / mminfo.mu;
   pd.nnu = mminfo.nb / mminfo.nu;
   pd.beta = beta;
   if (TA == AtlasNoTrans)
   {
      pd.incAm = pd.mb;
      pd.incAk = pd.kb * lda;
   }
   else
   {
      pd.incAm = pd.mb * lda;
      pd.incAk = pd.kb;
   }
   if (TB == AtlasNoTrans)
   {
      pd.incBk = pd.kb;
      pd.incBn = pd.nb * ldb;
   }
   else
   {
      pd.incBk = pd.kb * ldb;
      pd.incBn = pd.nb;
   }
   pd.incCn = nb*ldc;
   pd.incCm = mb;
   pd.KB0 = KB0;
/*
 * Allocate workspace, for now get double amount space required for first
 * recursively divided dimension!
 */
   i = (mr) ? nmblks+1 : nmblks;
   j = (nr) ? nnblks+1 : nnblks;
   k = (kr) ? nkblks+1 : nkblks;
   szA = i*k*pd.mbkb;
   szB = k*j*pd.kbnb;
   szC = i*j*pd.mbnb;
   vp = malloc(ATL_MulBySize(szA+mu + szB+nu + szC+mu*nu) + 3*ATL_Cachelen);
   ATL_assert(vp);
   a = ATL_AlignPtr(vp);
   b = a + szA;
   b = ATL_AlignPtr(b);
   c = b + szB;
   c = ATL_AlignPtr(c);
   flag = (v_cpA | v_cpB | v_cpC | v_nmA | v_nmB | v_nmC | v_lwC);
   if (!(mr|nr|kr))
      ammmRECf(&pd, nmblks, nnblks, nkblks, flag, A, B, C, a, b, c);
   else
   {
      flag |= (mr) ? v_mp : 0;
      flag |= (nr) ? v_np : 0;
      flag |= (kr) ? v_kp : 0;
      ammmREC(&pd, nmblks, nnblks, nkblks, flag, A, B, C, a, b, c);
   }
   free(vp);
}
