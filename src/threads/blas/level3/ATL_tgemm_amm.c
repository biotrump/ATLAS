#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"

/*
 * recurs on any standard GEMM interface routine, passed as a function
 * pointer in amm.
 */
int Mjoin(PATL,ammm_REC)
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
   ATL_CINT ldc,
   int (*amm)(enum ATLAS_TRANS,enum ATLAS_TRANS, ATL_CINT, ATL_CINT, ATL_CINT,
              const SCALAR, const TYPE*, ATL_CINT,  const TYPE*, ATL_CINT,
              const SCALAR, TYPE*, ATL_CINT)
)
{
   if (amm(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc))
   {
/*
 *    Stopping criteria in case something is horribly wrong
 */
      if (K <= ATL_AMM_MAXKB && M <= ATL_AMM_MAXMB && N <= ATL_AMM_MAXNB)
         return(Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                 beta, C, ldc));
/*
 *    Divide K first: it cuts space from both A & B
 */
      if (K+K >= Mmax(M,N))
      {
         const int KL = K>>1, KR = K-KL;
         ATL_assert(!Mjoin(PATL,ammm_REC)(TA, TB, M, N, KL, alpha, A, lda,
                                          B, ldb, beta, C, ldc, amm));
         A += (TA == AtlasNoTrans) ? KL*lda : KL;
         B += (TB == AtlasNoTrans) ? KL : KL*ldb;
         return(Mjoin(PATL,ammm_REC)(TA, TB, M, N, KR, alpha, A, lda,
                                     B, ldb, ATL_rone, C, ldc, amm));
      }
/*
 *    If M largest dim (twice K), cut it instead
 */
      else if (M >= N)
      {
         const int ML = M>>1, MR = M-ML;
         ATL_assert(!Mjoin(PATL,ammm_REC)(TA, TB, ML, N, K, alpha, A, lda,
                                          B, ldb, beta, C, ldc, amm));
         A += (TA == AtlasNoTrans) ? ML : ML*lda;
         return(Mjoin(PATL,ammm_REC)(TA, TB, MR, N, K, alpha, A, lda,
                                     B, ldb, beta, C+ML, ldc, amm));
      }
/*
 *    Otherwise, cut N
 */
      else
      {
         const int NL = N>>1, NR = N-NL;
         ATL_assert(!Mjoin(PATL,ammm_REC)(TA, TB, M, NL, K, alpha, A, lda,
                                          B, ldb, beta, C, ldc, amm));
         B += (TB == AtlasNoTrans) ? NL*ldb : NL;
         return(Mjoin(PATL,ammm_REC)(TA, TB, M, NR, K, alpha, A, lda,
                                     B, ldb, beta, C+NL*ldc, ldc, amm));
      }
   }
   return(0);
}

int Mjoin(PATL,tammm)
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
   if (N < ATL_AMM_MAXNB && K <= ATL_AMM_MAXKB && K > 2)
      return(Mjoin(PATL,tammm_tNK)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                   beta, C, ldc));
   else if (M <= ATL_AMM_MAXMB && N <= ATL_AMM_MAXNB)
      return(Mjoin(PATL,tammm_tMN)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                   beta, C, ldc));
   else if (M > ATL_AMM_MAXMB && K >= 16 && N > ATL_AMM_MAXNB)
      return(Mjoin(PATL,ammm_REC)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                  beta, C, ldc, Mjoin(PATL,tammm_G)));
   return(1);
}
