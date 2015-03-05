#include "atlas_misc.h"
#include "atlas_tlvl3.h"
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_use.h))

#if ATL_USE_AMM
/*
 * This function provides an estimate on max number of threads to use to
 * perform a access-major GEMM.
 */
size_t Mjoin(PATL,GetAmmmNthr)(ATL_CSZT M, ATL_CSZT N, ATL_CSZT K)
{
   size_t nnblks, nmblks, nkblks, p;

/*
 * On the XeonPHI, threads take a huge time to start up, so don't try
 * threading until we have a big problem.  We can hopefully fix this
 * later by changing to a thread pool model on the PHI.
 */
   #ifdef ATL_ARCH_XeonPHI
      if ((M*1e-6)*N*K < 8.0)
         return(1);
   #endif

   nmblks = (M >= ATL_AMM_66MB) ? M/ATL_AMM_66MB : 1;
   nnblks = (N >= ATL_AMM_66NB) ? N/ATL_AMM_66NB : 1;
   nkblks = (K >= ATL_AMM_66KB) ? K/ATL_AMM_66KB : 1;
/*
 * Any shape with two degenerate dimensions causes a lot of bus traffic,
 * with very little computation to overcome threading overheads,
 * so demand at least 32 blocks before parallelizing
 */
   if ((nmblks==1 && nnblks==1) || (nmblks==1 && nkblks==1) ||
       (nnblks==1 && nkblks==1))
      return((nnblks*nmblks*nkblks)>>5);
/*
 * If it is a rank-K update, ask to have 4 big blocks of C
 */
   if (K <= ATL_MAXM_RKK)
   {
      nnblks=N/ATL_MAXN_RKK, nmblks=M/ATL_MAXM_RKK;
      return((nnblks*nmblks)>>2);
   }

/*
 * By default, give everyone 32 blocks to compute; for square problems,
 * the number of blocks is cubic, so this should not meaningfully restrict
 * parallelism.
 */
   return((nmblks*nnblks*nkblks)>>5);
}
#endif

/*
 * ====================================================================
 * This function will eventually be generated, but for now just written
 * ====================================================================
 */
int Mjoin(PATL,tNumGemmThreads)(ATL_CINT M, ATL_CINT N, ATL_CINT K)
/*
 * RETURNS : estimate of how many threads will be used to run the problem,
 *           assuming we will actually do threading (i.e. THRESH is exceeded)
 *           0 is returned if this number is 1 or less.
 */
{
   #if ATL_USE_AMM
      size_t np;
      np = Mjoin(PATL,GetAmmmNthr)(M, N, K);
      return(np >= ATL_NTHREADS ? ATL_NTHREADS : np);
   #else
      double flops;
      int np;
      if (M < 4 || N < 4 || K < 4)
         return(0);
      flops = ((2.0*M)*N)*K;
      np = flops / ATL_TGEMM_PERTHR_MF;
      np = (np <= 1) ? 0 : np;
      return(np >= ATL_NTHREADS ? ATL_NTHREADS : np);
   #endif
}

int Mjoin(PATL,GemmWillThread)(ATL_CINT M, ATL_CINT N, ATL_CINT K)
/*
 * Returns: 0 if threshold FLOPS not achieved, rough # of threads used else
 */
{
   #if ATL_USE_AMM
      size_t np;
      np = Mjoin(PATL,GetAmmmNthr)(M, N, K);
      if (np < 2)
         return(0);
      return(Mmin(np, ATL_NTHREADS));
   #else
      if (((2.0*M)*N)*K < ATL_TGEMM_THRESH_MF)
         return(0);
      return(Mjoin(PATL,tNumGemmThreads)(M,N,K));
   #endif
}
