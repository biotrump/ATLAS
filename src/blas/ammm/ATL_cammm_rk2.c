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
   #ifdef DCPLX
      const int MB=512, NB=16;
   #else /* SCPLX */
      const int MB=512, NB=32;
   #endif
   void *vp;
   TYPE *x, *y, *w, *z;
   size_t j;
   const TYPE ONE[2] = {ATL_rone, ATL_rzero};
   ATL_CSZT lda2=lda+lda, ldb2=ldb+ldb, ldc2=ldc+ldc, incC = NB*ldc2;

/*
 * If beta is one, can handle by one call to ger2
 */
   if (SCALAR_IS_ONE(beta) && 0)
   {
      if (TA == AtlasNoTrans)
      {
         if (TB == AtlasNoTrans)
            Mjoin(PATL,ger2u)(M, N, alpha, A, 1, B, ldb, alpha, A+lda2, 1,
                              B+2, ldb, C, ldc);
         else if (TB == AtlasConj)
            Mjoin(PATL,ger2c)(M, N, alpha, A, 1, B, ldb, alpha, A+lda2, 1,
                              B+2, ldb, C, ldc);
         else if (TB == AtlasTrans)
            Mjoin(PATL,ger2u)(M, N, alpha, A, 1, B, 1, alpha, A+lda2, 1,
                              B+ldb2, 1, C, ldc);
         else /* if (TB == AtlasConjTrans) */
            Mjoin(PATL,ger2c)(M, N, alpha, A, 1, B, 1, alpha, A+lda2, 1,
                              B+ldb2, 1, C, ldc);
      }
      else if (TA == AtlasTrans)
      {
         if (TB == AtlasNoTrans)
            Mjoin(PATL,ger2u)(M, N, alpha, A, lda, B, ldb, alpha, A+2, lda,
                              B+2, ldb, C, ldc);
         else if (TB == AtlasConj)
            Mjoin(PATL,ger2c)(M, N, alpha, A, lda, B, ldb, alpha, A+2, lda,
                              B+2, ldb, C, ldc);
         else if (TB == AtlasTrans)
            Mjoin(PATL,ger2u)(M, N, alpha, A, lda, B, 1, alpha, A+2, lda,
                              B+ldb2, 1, C, ldc);
         else if (TB == AtlasConjTrans)
            Mjoin(PATL,ger2c)(M, N, alpha, A, lda, B, 1, alpha, A+2, lda,
                              B+ldb2, 1, C, ldc);
      }
/*
 *    If A must be conjugated, copy it
 */
      else  /* TA == AtlasConj || TA == AtlasConjTrans */
      {
         vp = malloc((ATL_MulBySize(M)+ATL_Cachelen)<<1);
         if (!vp)
           return(2);
         x = ATL_AlignPtr(vp);
         w = x + M + M;
         w = ATL_AlignPtr(w);
         if (SCALAR_IS_ONE(alpha))
         {
            if (TA == AtlasConj)
            {
               Mjoin(PATL,copyConj)(M, A, 1, x, 1);
               Mjoin(PATL,copyConj)(M, A+lda2, 1, w, 1);
            }
            else /* if (TA == AtlasConjTrans) */
            {
               Mjoin(PATL,copyConj)(M, A, lda, x, 1);
               Mjoin(PATL,copyConj)(M, A+2, lda, w, 1);
            }
         }
         else
         {
            if (TA == AtlasConj)
            {
               Mjoin(PATL,moveConj)(M, alpha, A, 1, x, 1);
               Mjoin(PATL,moveConj)(M, alpha, A+lda2, 1, w, 1);
            }
            else /* if (TA == AtlasConjTrans) */
            {
               Mjoin(PATL,moveConj)(M, alpha, A, lda, x, 1);
               Mjoin(PATL,moveConj)(M, alpha, A+2, lda, w, 1);
            }
         }
         if (TB == AtlasNoTrans)
            Mjoin(PATL,ger2u)(M, N, ONE, x, 1, B, ldb, ONE, w, 1,
                              B+2, ldb, C, ldc);
         else if (TB == AtlasConj)
            Mjoin(PATL,ger2c)(M, N, ONE, x, 1, B, ldb, ONE, w, 1,
                              B+2, ldb, C, ldc);
         else if (TB == AtlasTrans)
            Mjoin(PATL,ger2u)(M, N, ONE, x, 1, B, 1, ONE, w, 1,
                              B+ldb2, 1, C, ldc);
         else /* if (TB == AtlasConjTrans) */
            Mjoin(PATL,ger2c)(M, N, ONE, x, 1, B, 1, ONE, w, 1,
                              B+ldb2, 1, C, ldc);
         free(vp);
      }
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
   y = x + M + M;
   y = ATL_AlignPtr(y);
   w = y + N + N;
   w = ATL_AlignPtr(w);
   z = w + M + M;
   z = ATL_AlignPtr(z);
   if (TA == AtlasNoTrans)
   {
      Mjoin(PATL,copy)(M, A, 1, x, 1);
      Mjoin(PATL,copy)(M, A+lda2, 1, w, 1);
   }
   else if (TA == AtlasConj)
   {
      Mjoin(PATL,copyConj)(M, A, 1, x, 1);
      Mjoin(PATL,copyConj)(M, A+lda2, 1, w, 1);
   }
   else if (TA == AtlasTrans)
   {
      Mjoin(PATL,copy)(M, A, lda, x, 1);
      Mjoin(PATL,copy)(M, A+2, lda, w, 1);
   }
   else if (TA == AtlasConjTrans)
   {
      Mjoin(PATL,copyConj)(M, A, lda, x, 1);
      Mjoin(PATL,copyConj)(M, A+2, lda, w, 1);
   }
   if (SCALAR_IS_ONE(alpha))
   {
      if (TB == AtlasNoTrans)
      {
         Mjoin(PATL,copy)(N, B, ldb, y, 1);
         Mjoin(PATL,copy)(N, B+2, ldb, z, 1);
      }
      else if (TB == AtlasConj)
      {
         Mjoin(PATL,copyConj)(N, B, ldb, y, 1);
         Mjoin(PATL,copyConj)(N, B+2, ldb, z, 1);
      }
      else if (TB == AtlasTrans)
      {
         Mjoin(PATL,copy)(N, B, 1, y, 1);
         Mjoin(PATL,copy)(N, B+ldb2, 1, z, 1);
      }
      else if (TB == AtlasConjTrans)
      {
         Mjoin(PATL,copyConj)(N, B, 1, y, 1);
         Mjoin(PATL,copyConj)(N, B+ldb2, 1, z, 1);
      }
   }
   else  /* alpha non-one; must apply */
   {
      if (TB == AtlasNoTrans)
      {
         Mjoin(PATL,cpsc)(N, alpha, B, ldb, y, 1);
         Mjoin(PATL,cpsc)(N, alpha, B+2, ldb, z, 1);
      }
      else if (TB == AtlasConj)
      {
         Mjoin(PATL,moveConj)(N, alpha, B, ldb, y, 1);
         Mjoin(PATL,moveConj)(N, alpha, B+2, ldb, z, 1);
      }
      else if (TB == AtlasTrans)
      {
         Mjoin(PATL,cpsc)(N, alpha, B, 1, y, 1);
         Mjoin(PATL,cpsc)(N, alpha, B+ldb2, 1, z, 1);
      }
      else /* if (TB == AtlasConjTrans) */
      {
         Mjoin(PATL,moveConj)(N, alpha, B, 1, y, 1);
         Mjoin(PATL,moveConj)(N, alpha, B+ldb2, 1, z, 1);
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
         Mjoin(PATL,gescal)(mb, nb, beta, C+i+i, ldc);
         Mjoin(PATL,ger2u)(mb, nb, ONE, x+i+i, 1, y+j+j, 1, ONE,
                           w+i+i, 1, z+j+j, 1, C+i+i, ldc);
      }
   }
   free(vp);
   return(0);
}
