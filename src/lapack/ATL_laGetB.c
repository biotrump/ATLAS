#include "atlas_misc.h"
#include Mstr(Mjoin(AMM_UPR,_blk.h))
#include Mstr(Mjoin(AMM_PRE,_sum.h))
#include "atlas_lapack.h"
/*
 * Take a wild-ass guess as to a good blocking factor.
 * M, N, K dim GEMM called wt, 0 means that dim is blocked.
 * Type of operation is given in the HINT:
 *    0 : panel op is fast (eg. L3-BLAS based) w/o extra flops (eg. LU)
 *    1 : panel op is slow (eg, L1/2-BLAS based) without extra flops
 *    2 : panel op is fast, with extra flops (eg. QR)
 *    4 : panel op is slow (eg, L1/2-BLAS based) with extra flops
 * Right now, this guys assumes a rank-K update shape regardless of input;
 * it should be replaced with something smarter once we've finalized the
 * amm blocking strategy.  All the computation shouldn't disguise the fact
 * there is no real intelligence here.
 */
#ifdef ATL_CAMM_MAXINDX
   #define MAXIDX ATL_CAMM_MAXINDX-1
#else
   #define MAXIDX ATL_AMM_NCASES-1
#endif
int Mjoin(PATL,laGetB)
(
   int M, /* 0 if dim is blocked, else max size of rows of C in GEMM call */
   int N, /* 0 if dim is blocked, else max size of cols of C in GEMM call */
   int K, /* 0 if dim is blked, else max size of common A/B dim in GEMM call */
   int HINT   /* type of operation to be blocked */
)
{
   #ifdef DCPLX
      const int maxB = (HINT < 2) ? ATL_AMM_KBs[MAXIDX] : 44;
   #elif defined(SCPLX) || defined(DREAL)
      const int maxB = (HINT < 2) ? ATL_AMM_KBs[MAXIDX] : 60;
   #else /* SREAL */
      const int maxB = (HINT < 2) ? ATL_AMM_KBs[MAXIDX] : 80;
   #endif
   int minblks, majSH, minSH=0, D = (M) ? M : N;
   int B, i;

   D = (D) ? D : K;
   B = ATL_AMM_KBs[0];
   if ((B<<majSH)+(B<<minSH) >= D)  /* quick exit for tiny problem */
   {
      if ((B<<2) < D)
         return(1);
      return(B);
   }
   if (HINT == 4)       /* slow, extra flops */
   {
      minblks = 48;    /* (32+16) */
      majSH = 5;
      minSH = 4;
   }
   else if (HINT == 1)  /* slow, no extra flops */
   {
      minblks = 24;    /* 16+8 */
      majSH = 4;
      minSH = 3;
   }
   else if (HINT == 2)  /* fast, extra flops */
   {
      minblks = 12;
      majSH = 3;
      minSH = 2;
   }
   else                 /* fast, no extra flops */
   {
      minblks = 9;
      majSH = 3;
      minSH = 1;
   }
   for (i=MAXIDX; i > 0; i--)
      if (B <= maxB)
         break;
   if ((B<<majSH)+(B<<minSH) <= D)  /* large problem, can use maxB */
      return(B);

   for (i--; i > 0; i--)
   {
      B = ATL_AMM_KBs[i];
      if ((B<<majSH)+(B<<minSH) <= D)
         return(B);
   }
   return(ATL_AMM_KBs[0]);
}
