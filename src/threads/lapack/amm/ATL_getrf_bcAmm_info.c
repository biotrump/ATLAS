#include "atlas_misc.h"
#include "atlas_pthreads.h"
#include "atlas_amm.h"
#define rpref Mjoin(Mjoin(Mjoin(atlas_,PRE),u0),amm)
#define cpref Mjoin(Mjoin(Mjoin(atlas_,PRE),u0),amm)
#include Mstr(Mjoin(rpref,_perf.h))
#include Mstr(Mjoin(rpref,_blk.h))

#if !defined(ATL_NTHREADS) || !defined(ATL_BC_RMAX)
   #ifndef ATL_NTHREADS
      #error "ATL_NTHREADS must be defined!"
   #elif ATL_NTHREADS < 8
      #define ATL_BC_RMAX 1
   #elif ATL_NTHREADS < 18
      #define ATL_BC_RMAX 2
   #elif ATL_NTHREADS < 32
      #define ATL_BC_RMAX 3
   #elif ATL_NTHREADS < 50
      #define ATL_BC_RMAX 4
   #elif ATL_NTHREADS < 72
      #define ATL_BC_RMAX 5
   #elif ATL_NTHREADS < 98
      #define ATL_BC_RMAX 6
   #elif ATL_NTHREADS < 128
      #define ATL_BC_RMAX 7
   #elif ATL_NTHREADS < 162
      #define ATL_BC_RMAX 8
   #elif ATL_NTHREADS < 200
      #define ATL_BC_RMAX 9
   #elif ATL_NTHREADS < 242
      #define ATL_BC_RMAX 10
   #elif ATL_NTHREADS < 288
      #define ATL_BC_RMAX 11
   #elif ATL_NTHREADS < 338
      #define ATL_BC_RMAX 12
   #elif ATL_NTHREADS < 392
      #define ATL_BC_RMAX 13
   #elif ATL_NTHREADS < 450
      #define ATL_BC_RMAX 14
   #elif ATL_NTHREADS < 512
      #define ATL_BC_RMAX 15
   #elif ATL_NTHREADS < 578
      #define ATL_BC_RMAX 16
   #elif ATL_NTHREADS < 648
      #define ATL_BC_RMAX 17
   #elif ATL_NTHREADS < 722
      #define ATL_BC_RMAX 18
   #elif ATL_NTHREADS < 800
      #define ATL_BC_RMAX 19
   #elif ATL_NTHREADS < 882
      #define ATL_BC_RMAX 20
   #elif ATL_NTHREADS < 968
      #define ATL_BC_RMAX 21
   #elif ATL_NTHREADS < 1058
      #define ATL_BC_RMAX 22
   #endif
#endif
static INLINE int findmaxR(int p)
{
   int r;
   for (p >>= 1, r=1; r*r < p; r++);
   if (r*r != p)  /* if not equal */
      r--;        /* round down */
   return(r);
}
#ifdef TCPLX
   #define GetTRSMflops(fc_, n_, r_) \
      fc_ += (4.0*(r_))*(n_) * ((n_) - 1.0)
   #define GetGEMMflops(fc_, m_, n_, k_) \
      fc_ += ((8.0*m_)*(n_))*(k_);
#else
   #define GetTRSMflops(fc_, n_, r_) \
      fc_ += ((double)(r_))*(n_) * ((n_) - 1.0)
   #define GetGEMMflops(fc_, m_, n_, k_) \
      fc_ += ((2.0*m_)*(n_))*(k_);
#endif
static INLINE double GetLUflops(int M, int N)
{
   double m=M, n=N, k = (M <= N)?m:n, ops;
   ops = m*n - 0.5*(m+n)*(k+1) + (1.0/6.0)*(k+1)*(k+k+1);
   ops = k * ( ops+ops + m - 0.5*(k+1) );
   #ifdef TCPLX
      ops *= 4.0;
   #endif
   return(ops);
}


int SumL3Flops(ATL_CUINT M, ATL_CUINT N, ATL_CUINT NB, double *MFL, double *SFL)
{
   double mfl, sfl;
   ATL_CUINT MN = (N <= M) ? N : M;
   ATL_UINT j;
   int np;

   mfl = sfl = 0.0;
   np = 1;
   for (j=0; j < MN; j += NB, np++)
   {
      int jb = MN-j, n, m;
      jb = (jb > NB) ? NB : jb;

      n = N-j-jb;
      if (n > 0)
      {
         GetTRSMflops(sfl, jb, n);
         m = M-j-jb;
         if (m > 0)
            GetGEMMflops(mfl, m, n, jb);
      }
   }
   *MFL = mfl;
   *SFL = sfl;
   return(np);
}
#define MINROWBLKS 128
int Mjoin(PATL,getrf_bcAmm_info)   /* RETURNS: suggested amm index */
(
   ATL_INT M,  /* rows of matrix */
   ATL_INT N,  /* columns of matrix */
   int *NB,    /* blocking to use */
   int *R,     /* OUTPUT: suggested R of RxC process grid */
   int *C
)
{
   int idx=0, P, i, NNB, r, c, nb;
   double tfl, pfl, mfl, rfl;   /* total, panel, matmul, trsm */
   #ifndef TRSM_OVER_GEMM
      double trsmM = 0.25;  /* say nb-trsm runs a quarter gemm's speed */
   #else
      double trsmM = TRSM_OVER_GEMM;
   #endif
   #ifndef GETF2_OVER_GEMM
      double panM = 0.25;  /* say panal fact runs quarter gemm's speed */
   #else
      double panM = GETF2_OVER_GEMM;
   #endif
   int maxR;

   P = Mjoin(PATL,tNumGemmThreads)(M, N>>1, N);
   if (P < 2)
   {
      *R = 1;
      *C = 1;
      *NB = 4;
      return(0);
   }
   if (P < ATL_NTHREADS)
      maxR = findmaxR(P);
   #ifndef ATL_BC_RMAX
      maxR = findmaxR(ATL_NTHREADS);
   #else
      maxR = ATL_BC_RMAX;
   #endif
   i = (N>>7);
   maxR = (maxR <= i) ? maxR : i;
   while ((P/maxR)*maxR != P)
      maxR--;
   r = *R = maxR;
   c = *C = P / maxR;
/*
 * Replace this with a model; it presently sucks
 */
   for (i=ATL_UAMM_NCASES-1; i; i--)
   {
      nb = ATL_UAMM_NBs[i];
      if (N > (c<<1)*nb)
         break;
   }
   idx = i;
   *NB = nb;
/*
 * For each panel factorization, think we should again compute % of gemm
 * flops, where we say NB = N / log2(N).  Now, we predict the GEMM code
 * will run at the speed a GEMM of that NB * R (time /R).  Predict the
 * rest of the panel gets no speedup.
 * Predict TRSM will gain parallelism only from C, but will be perfect.
 * Predict GEMM gets speedup wt R*C, but that row has slowdown (copy).
 */
#if 0
   tfl = GetLUflops(M, N);
   pfl = SumPanelFlops(M, N, nb, &NNB);
#endif
   printf("IDXAMM=%d: nb=%d, R=%d C=%d\n", idx, *NB, *R, *C);
   return(idx);
}
