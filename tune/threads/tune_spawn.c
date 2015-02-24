#include "atlas_threads.h"
#include "atlas_misc.h"


void TuneDoWork(ATL_LAUNCHSTRUCT_t *lp, void *vp)
/*
 * Use volatile array to check in, and then quit (cache-speed barrier)
 */
{
   ATL_thread_t *tp = vp;
   const int nthr = tp->P, rank=tp->rank;
   int i;
   volatile int *donearr=lp->opstruct;

   donearr[rank] = 1;
   for (i=0; i < nthr; i++)
      while(!donearr[i]);
}

void PrintUsage(char *exe)
{
   fprintf(stderr, "USAGE: %s [-r <reps>] -W [which]\n", exe);
   fprintf(stderr,
   "   which: bitfield, 1st bit is dyn, 2nd is lg2, 3rd is linear\n");
   exit(-1);
}

int GetFlags(int nargs, char **args, int *which)
{
   int i, reps=1;

   *which = 7;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'r':
         if (++i >= nargs)
            PrintUsage(args[0]);
         reps = atoi(args[i]);
         break;
      case 'W':
         if (++i >= nargs)
            PrintUsage(args[0]);
         *which = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
   return(reps);
}
int main(int nargs, char **args)
{
   int i, k, nreps = 200, which;
   double t0, tlin, tlg2, tdyn, trnk;
   volatile int done[ATL_NTHREADS];

   tlg2 = tdyn = tlin = 0.0;
   nreps = GetFlags(nargs, args, &which);


   printf("FINDING SPEED OF CREATE/BARRIER/JOIN USING %d REPITITIONS:\n",
          nreps);
   if (which & 1)
   {
      t0 = ATL_walltime();
      for (k=0; k < nreps; k++)
      {
         for (i=0; i < ATL_NTHREADS; i++) done[i] = 0;
         ATL_goparallel_dyn(ATL_NTHREADS, TuneDoWork, done, NULL);
      }
      tdyn = ATL_walltime() - t0;
      printf("   dyn time = %e\n", (float)tdyn);
   }

   if (which & 2)
   {
      t0 = ATL_walltime();
      for (k=0; k < nreps; k++)
      {
         for (i=0; i < ATL_NTHREADS; i++) done[i] = 0;
         ATL_goparallel_log2(ATL_NTHREADS, TuneDoWork, done, NULL);
      }
      tlg2 = ATL_walltime() - t0;
      printf("   lg2 time = %e\n", (float)tlg2);
   }

   if (which & 4)
   {
      t0 = ATL_walltime();
      for (k=0; k < nreps; k++)
      {
         for (i=0; i < ATL_NTHREADS; i++) done[i] = 0;
         ATL_goparallel_lin(ATL_NTHREADS, TuneDoWork, done, NULL);
      }
      tlin = ATL_walltime() - t0;
      printf("   lin time = %e\n", (float)tlin);
   }
   #if 0
   if (which & 8)
   {
      t0 = ATL_walltime();
      for (k=0; k < nreps; k++)
      {
         for (i=0; i < ATL_NTHREADS; i++) done[i] = 0;
         ATL_goparallel_prank(ATL_NTHREADS, TuneDoWork_gp, ta, NULL);
      }
      trnk = ATL_walltime() - t0;
      printf("   rnk time = %e\n", (float)trnk);
   }
   #endif
   if ((which | 7) == which)
      printf("DYNAMIC is %.2f%% of LINEAR and %.2f%% of LOG2 SPEED.\n",
             (tdyn/tlin)*100.0, (tdyn/tlg2)*100.0);
   if ((which & 1) && (which & 8))
      printf("rank dynamic is %.2f%% of affinity dynamic\n", (trnk/tdyn)*100.0);
   return(0);
}
