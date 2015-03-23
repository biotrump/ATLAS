#include "atlas_misc.h"
#include "atlas_lvl2.h"
#include "atlas_level1.h"
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_sum.h))
static int ATL_ammm
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
/*
 * Just do a scale and return
 */
   if (SCALAR_IS_ZERO(alpha) || !K)
   {
      if (SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gezero)(M, N, C, ldc);
      else if (!SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);
      return(0);
   }
/*
 * Scope for degenerate cases that should call Level-2 BLAS; these
 * routines assert they work, since their workspace is O(N) and so are
 * not allowed to fail.
 */
   if (K == 1)  /* really a GER */
   {
      if (!SCALAR_IS_ONE(beta))
      {
         int i;
         ATL_CSZT incA = ((TA == AtlasNoTrans) ? 1 : lda);
         ATL_CSZT incB = ((TB == AtlasNoTrans) ? ldb : 1);
         for (i=0; i < N; i++, B += incB, C += ldc)
            Mjoin(PATL,axpby)(M, alpha * *B,  A, incA, beta, C, 1);
      }
      else
         Mjoin(PATL,ger)(M, N, alpha, A, (TA == AtlasNoTrans) ? 1 : lda,
                         B, (TB == AtlasNoTrans) ? ldb : 1, C, ldc);
      return(0);
   }
   if (K == 2)
      return(Mjoin(PATL,ammm_rk2)(TA, TB, M, N, alpha, A, lda, B, ldb,
                                  beta, C, ldc));
   if (N == 1)  /* Really GEMV with A as matrix, A & C as vectors */
   {
      if (TA == AtlasNoTrans)
         Mjoin(PATL,gemv)(AtlasNoTrans, M, K, alpha, A, lda, B,
                          (TB == AtlasNoTrans) ? 1:ldb, beta, C, 1);
      else
         Mjoin(PATL,gemv)(AtlasTrans, K, M, alpha, A, lda, B,
                          (TB == AtlasNoTrans) ? 1:ldb, beta, C, 1);
      return(0);
   }
   if (M == 1)  /* Really GEMV with B as matrix, A & C as vectors */
   {
      if (TB == AtlasNoTrans)
         Mjoin(PATL,gemv)(AtlasTrans, K, N, alpha, B, ldb, A,
                          (TA == AtlasNoTrans) ? lda:1, beta, C, ldc);
      else
         Mjoin(PATL,gemv)(AtlasNoTrans, N, K, alpha, B, ldb, A,
                          (TA == AtlasNoTrans) ? lda:1, beta, C, ldc);
      return(0);
   }
/*
 * Special case mainly for LU, where K==4, N==4 beta=1.0, TA==AtlasNoTrans;
 * Can do a no-copy update with only one loop.
 */
   if (K == 4 && N == 4 && TA==AtlasNoTrans && SCALAR_IS_ONE(beta))
   {
      int Mjoin(PATL,rk4n4)(enum ATLAS_TRANS,ATL_CSZT,const SCALAR,
          const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,TYPE*,ATL_CSZT);
      if (!Mjoin(PATL,rk4n4)(TB, M, alpha, A, lda, B, ldb, C, ldc))
         return(0);
   }
/*
 * 1-block special case code can return w/o doing op if it thinks
 * rank-K would be faster
 */
   if (M <= ATL_AMM_MAXMB && N <= ATL_AMM_MAXNB && K <= ATL_AMM_MAXKB)
      if (!Mjoin(PATL,ammm_1b)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                               beta, C, ldc))
         return(0);
/*
 * Rank-K can fail to allocate space, so return success/failure
 */
   if (K > 2 && K <= ATL_MAXK_RKK)
      return(Mjoin(PATL,ammm_rkK)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                  beta, C, ldc));
/*
 * Handle case that is really an inner product shape (M<=MB, N<=NB, large K)
 * This case not allowed to fail since it requires only 3*NB^2 workspace.
 */
   if (M <= ATL_AMM_MAXMB && N <= ATL_AMM_MAXNB)
      return(Mjoin(PATL,ammm_IP)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                 beta, C, ldc));
/*
 * If B/C have only one column panel, call special low-workspace (3NB^3)
 * code for additional performance.  This shape occurs in left-looking LU.
 */
   if (N <= ATL_AMM_MAXNB)
      return(Mjoin(PATL,ammm_tN)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                 beta, C, ldc));
/*
 * Next two loop orderings are general case, so use whichever uses least
 * workspace
 */
   if (M >= N)
      return(Mjoin(PATL,ammmMNK)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                 beta, C, ldc));
/*
 * This guy tries to allocate (M+NB)*K + NB^2 worskpace, so recursion
 * may be needed to keep it within allotted memory.
 */
   return(Mjoin(PATL,ammmNMK)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                              beta, C, ldc));
}
/*
 * Recur to get K below this value; this puts a ceiling on workspace and
 * usually improves performance (in huge problems, reduces TLB pressure)
 */
#define ATL_MAX_RK 3000


/*
 * This routine uses recursion to cut the dimensions of the matrices until
 * workspace requirements are low enough that a call to ATL_ammm succeeds
 */
void Mjoin(PATL,ammm)
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
   if (!M || !N)
      return;
/*
 * Cases where all we must do is possibly scale and return
 */
   if (SCALAR_IS_ZERO(alpha) || !K)
   {
      if (SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gezero)(M, N, C, ldc);
      else if (!SCALAR_IS_ONE(beta))
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);
      return;
   }
/*
 * Our stopping criteria is if ATL_ammm signals success in mallocing mem
 */
   if (K <= ATL_MAX_RK)
   {
      if(!ATL_ammm(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc))
         return;
   }
/*
 * =========================================================================
 * Otherwise, problem too large, so we'll recursively divide its largest dim
 * =========================================================================
 */
/*
 * if K is tied for largest, cut it, since it reduces size of A & B
 * NOTE: C always uses only NB^2 workspace, so only A/B matters.
 */
   if (K > ATL_MAX_RK || (K >= N && K >= M))
   {
      const size_t kL=(K>>4)<<3, kR=K-kL;

      Mjoin(PATL,ammm)(TA, TB, M, N, kL, alpha, A, lda, B, ldb, beta, C, ldc);
      if (TA == AtlasNoTrans || TA == AtlasConj)
         A += (lda*kL)SHIFT;
      else
         A += kL SHIFT;
      if (TB == AtlasNoTrans)
         B += kL SHIFT;
      else
         B += (ldb*kL) SHIFT;
      Mjoin(PATL,ammm)(TA, TB, M, N, kR, alpha, A, lda, B, ldb,
                       ATL_rone, C, ldc);
   }
   else if (N >= M)  /* cutting N */
   {
      const size_t nL = (N>>1), nR = N-nL;
      Mjoin(PATL,ammm)(TA, TB, M, nL, K, alpha, A, lda, B, ldb, beta, C, ldc);
      if (TB == AtlasNoTrans || TB == AtlasConj)
         B += (ldb*nL)SHIFT;
      else
         B += nL SHIFT;
      C += (ldc*nL)SHIFT;
      Mjoin(PATL,ammm)(TA, TB, M, nR, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
   else  /* cutting M */
   {
      const size_t mL = (M>>1), mR = M-mL;
      Mjoin(PATL,ammm)(TA, TB, mL, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      if (TA == AtlasNoTrans || TA == AtlasConj)
         A += mL SHIFT;
      else
         A += (mL*lda)SHIFT;
      C += mL SHIFT;
      Mjoin(PATL,ammm)(TA, TB, mR, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
}

void Mjoin(PATL,gemm)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                      const int M, const int N, const int K, const SCALAR alpha,
                      const TYPE *A, const int lda, const TYPE *B,
                      const int ldb, const SCALAR beta, TYPE *C, const int ldc)
{
   Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
