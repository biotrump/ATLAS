#include "atlas_misc.h"
#include <assert.h>

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -n <spclen>: workspace to pass\n");
   fprintf(stderr, "   -I <its> : iterations to pass\n");
   fprintf(stderr, "   -# <ntimes> : set # of times to time kernel\n");
   exit(ierr ? ierr : -1);
}
size_t GetFlags(int nargs, char **args, size_t *N, int *NREP)
{
   size_t nits = 10000;
   int i;

   *N = 128;
   *NREP = 3;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'I':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         nits = atoll(args[i]);
         break;
      case 'n':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *N = atoll(args[i]);
         break;
      case '#':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *NREP = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   return(nits);
}

double ATL_walltime(void);
void RunKern(size_t nits, void *vp);

int main (int nargs, char **args)
{
   size_t n, nits, i, ifl;
   double t0, mf, mfB;
   void *vp;
   char *cp;
   size_t *lp;
   int nrep;

   nits = GetFlags(nargs, args, &n, &nrep);
   vp = malloc(n+32);
   cp = ATL_AlignPtr(vp);
   lp = (size_t*)cp;
   for (i=0; i < n; i++)
      cp[i] = 0;

   mfB=0;
   for (i=0; i < nrep; i++)
   {
      t0 = ATL_walltime();
      RunKern(nits, (void*)cp);
      t0 = ATL_walltime() - t0;
      ifl = lp[0];
      mf = (((double)ifl)*nits) / (t0 * 1.0e6);
      if (mf > mfB)
         mfB = mf;
      printf("MFLOPS=%.2f\n", mf);
   }
   printf("\nBEST = %.2f\n\n", mfB);

   free(vp);
   return(0);
}
