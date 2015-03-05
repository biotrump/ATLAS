#include "atlas_misc.h"
#include "atlas_lvl2.h"
#include "atlas_level1.h"
/*
 * This is special-case code that handles rank-2 update by calling GER2
 */
int Mjoin(PATL,ammm_rk2)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
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
   #ifdef DREAL
      const int MB=512, NB=32;
   #else /* SREAL */
      const int MB=512, NB=64;
   #endif
   void *vp;
   TYPE *x, *y, *w, *z;
   size_t j;
   ATL_CSZT incC = NB*ldc;

/*
 * If beta is one, can handle by one call to ger2
 */
   if (SCALAR_IS_ONE(beta))
   {
      if (TA == AtlasNoTrans)
      {
         if (TB == AtlasNoTrans)
            Mjoin(PATL,ger2)(M, N, alpha, A, 1, B, ldb, alpha, A+lda, 1,
                             B+1, ldb, C, ldc);
         else
            Mjoin(PATL,ger2)(M, N, alpha, A, 1, B, 1, alpha, A+lda, 1,
                             B+ldb, 1, C, ldc);
      }
      else if (TB == AtlasNoTrans)
         Mjoin(PATL,ger2)(M, N, alpha, A, lda, B, ldb, alpha, A+1, lda,
                          B+1, ldb, C, ldc);
      else
         Mjoin(PATL,ger2)(M, N, alpha, A, lda, B, 1, alpha, A+1, lda,
                          B+ldb, 1, C, ldc);
      return(0);
   }
/*
 * Later on, do smart think like copy only MB/NB at a time, and don't copy
 * at all if vectors are contiguous, but right now, always do copy up-front
 * so loop does not have to worry about TA/TB; this is a O(N) cost in N^2 alg
 */
   vp = malloc(2*ATL_MulBySize(M+N)+4*ATL_Cachelen);
   if (!vp)
      return(1);
   x = ATL_AlignPtr(vp);
   y = x + M;
   y = ATL_AlignPtr(y);
   w = y + N;
   w = ATL_AlignPtr(w);
   z = w + M;
   z = ATL_AlignPtr(z);
   if (TA == AtlasNoTrans)
   {
      Mjoin(PATL,copy)(M, A, 1, x, 1);
      Mjoin(PATL,copy)(M, A+lda, 1, w, 1);
   }
   else
   {
      Mjoin(PATL,copy)(M, A, lda, x, 1);
      Mjoin(PATL,copy)(M, A+1, lda, w, 1);
   }
   if (SCALAR_IS_ONE(alpha))
   {
      if (TB == AtlasNoTrans)
      {
         Mjoin(PATL,copy)(N, B, ldb, y, 1);
         Mjoin(PATL,copy)(N, B+1, ldb, z, 1);
      }
      else
      {
         Mjoin(PATL,copy)(N, B, 1, y, 1);
         Mjoin(PATL,copy)(N, B+ldb, 1, z, 1);
      }
   }
   else
   {
      if (TB == AtlasNoTrans)
      {
         Mjoin(PATL,cpsc)(N, alpha, B, ldb, y, 1);
         Mjoin(PATL,cpsc)(N, alpha, B+1, ldb, z, 1);
      }
      else
      {
         Mjoin(PATL,cpsc)(N, alpha, B, 1, y, 1);
         Mjoin(PATL,cpsc)(N, alpha, B+ldb, 1, z, 1);
      }
   }
   for (j=0; j < N; j += NB, C += incC)
   {
      size_t i, nb = N-j;
      nb = (nb >= NB) ? NB : nb;
      for (i=0; i < M; i += MB)
      {
         size_t mb = M-i;
         mb = (mb >= MB) ? MB : mb;
         Mjoin(PATL,gescal)(mb, nb, beta, C+i, ldc);
         Mjoin(PATL,ger2)(mb, nb, ATL_rone, x+i, 1, y+j, 1, ATL_rone,
                          w+i, 1, z+j, 1, C+i, ldc);
      }
   }
   free(vp);
   return(0);
}
