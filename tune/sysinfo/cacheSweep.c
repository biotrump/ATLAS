#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#if defined(CPUTIME)
   #define time00 ATL_cputime
#else
   #define time00 ATL_walltime
#endif
double time00();

void InitSpace(size_t N, double *dp)
{
   size_t i;
   for (i=0; i < N; i++)
      dp[i] = 0.0;
}

double ReadSpace(size_t N, double *dp)
{
   size_t i, N16 = (N>>4)<<4;
   double d=0.0, d1=0.0, d2=0.0, d3=0.0, d4=0.0, d5=0.0, d6=0.0, d7=0.0;
   for (i=0; i < N16; i += 16)
   {
      d += dp[i];
      d1 += dp[i+2];
      d2 += dp[i+4];
      d3 += dp[i+6];
      d4 += dp[i+8];
      d5 += dp[i+10];
      d6 += dp[i+12];
      d7 += dp[i+14];
   }
   return(d+d1+d2+d3);
}

double *DoCacheSearch(int fixreps, int timreps, size_t minsz, size_t maxsz,
                      int *N, size_t **SZS)
{
   size_t i, n=0, nrep, K;
   void *vp;
   double *dp, sum, *times, t0, t1;
   size_t *szs;
   printf("SIZE (KB)            TIME     BW (MB/s)\n");
   printf("=========  ==============  ============\n");
/*
 * Get space big enough to hold largest size
 */
   vp = malloc(maxsz + 512);
   assert(vp);
/*
 * Figure out how many timings we are going to do, and allocate time and
 * dimension arrays
 */
   for (i=maxsz; i >= minsz; i >>= 1)
      n++;
   *N = n;
   *SZS = szs = malloc(sizeof(size_t)*n);
   times = malloc(sizeof(double)*n);
   assert(szs && times);

   dp = (double*)((((size_t)vp)>>9)<<9);
   K = maxsz / sizeof(double);
   nrep = timreps;
   InitSpace(K, dp);

   for (n=0,i=maxsz; i >= minsz; i >>= 1, n++, nrep <<= 1)
   {
      int j;
      szs[n] = i;
      K = i / sizeof(double);
/*
 *    Get memory area allocated to the cache as much as probability allows
 *    by repeatedly sweeping it
 */
      sum = 0.0;
      for (j=0; j < fixreps; j++)
         sum += ReadSpace(K, dp);
      assert(sum == 0.0);   /* keeps compiler from eliminating sweeps */
      t0 = time00();
      for (j=0; j < nrep; j++)
         sum += ReadSpace(K, dp);
      times[n] = time00() - t0;
      assert(sum == 0.0);   /* keeps compiler from eliminating sweeps */
      t0 = ((i/1024)/1024.0)*nrep / times[n];
      printf("%9d  %14e  %12.0f\n", (int)(i / 1024), times[n], t0);
   }

   return(times);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -r <#> : # of reads of max size to perform\n");
   fprintf(stderr, "   -f <#> : set trip count to fix mem to cache\n");
   fprintf(stderr, "   -M <#> : max size to time (KB) \n");
   fprintf(stderr, "   -m <#> : min size to time (KB) \n");
   exit(ierr ? ierr : -1);
}

int GetFlags(int nargs, char **args, size_t *MINSZ, size_t *MAXSZ, int *FREP)
{
   int timrep=512, i;
   *MINSZ = 4;
   *MAXSZ = 8*1024;
   *FREP = 16;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'M':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *MAXSZ = atol(args[i]);
         break;
      case 'm':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *MINSZ = atol(args[i]);
         break;
      case 'f':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *FREP = atoi(args[i]);
         break;
      case 'r':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         timrep = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   return(timrep);
}
int main(int nargs, char **args)
{
   int fixrep, timrep, N;
   size_t minsz, maxsz, *szs;
   double *times;
   timrep = GetFlags(nargs, args, &minsz, &maxsz, &fixrep);
   minsz *= 1024;
   maxsz *= 1024;
   printf("MAXSZ=%lu, MINSZ=%lu, fixrep=%d, timrep=%d\n\n",
          (unsigned long)maxsz, (unsigned long)minsz, fixrep, timrep);
   times = DoCacheSearch(fixrep, timrep, minsz, maxsz, &N, &szs);
   free(times);
   free(szs);
   return(0);
}
