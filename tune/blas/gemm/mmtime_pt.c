#define _GNU_SOURCE 1 /* what manpage says you need to get CPU_SET */
#define __USE_GNU   1 /* what you actually have to set on linuxes I've seen */
#include <sched.h>    /* must include this before pthreads */
#include <pthread.h>
#include "atlas_ptmisc.h"
#include <assert.h>
#include <string.h>
#define dumb_rand() ( 0.5 - ((double)rand())/((double)RAND_MAX) )

#ifndef KMM
   #define KMM ATL_USERMM
#endif

#define CINT const int

#ifdef ATL_NEWTIME
   void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
            const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
            const TYPE *pCn);
#else
   void KMM(const int, const int, const int, const SCALAR, const TYPE*,
            const int, const TYPE*, const int, const SCALAR, TYPE*, const int);
#endif

struct kmm_struct{
   int mb, nb, kb;                      /* C: mbxnb, At: kbxmb, B: kbXnb */
   #ifdef ATL_NEWTIME
   int mu, nu, ku;                      /* needed to compute mblks/nblks */
   #endif
   int movA, movB, movC;                /* which mat move in flush array? */
   int FLSIZE;                          /* min area to move in in bytes */
   int reps;                            /* # calls to kmm in one timing */
   int LDC;                             /* what should ldc be set to? */
   int iam;                             /* thread rank */
   int p;                               /* total number of threads */
   int *pids;                           /* IDs of processors for affinity */
   double mf;                           /* mflop returned by timing */
};

double GetKmmMflop
(
   CINT mb, CINT nb, CINT kb,           /* C: mbxnb, At: kbxmb, B: kbXnb */
   #ifdef ATL_NEWTIME
      CINT mu, CINT nu, CINT ku,
   #endif
   CINT movA, CINT movB, CINT movC,     /* which mat move in flush array? */
   int FLSIZE,                          /* min area to move in in bytes */
   CINT reps,                           /* # calls to kmm in one timing */
   CINT LDC                             /* what should ldc be set to? */
)
/*
 * Returns MFLOP rate of matmul kernel KMM
 * LDC: if (LDC == 0), then set ldc=MB for timings.
 *      if (LDC != 0 && movC != 0), then ldc= col length in move space
 *      else ldc = LDC;
 *
 */
{
   #ifdef ATL_NEWTIME
      CINT mblks = mb/mu, nblks = nb/nu;
   #endif
   const int NOMOVE = !(movA|movB|movC);
   int ldc, setsz, nset, i, j, incA, incB, incC, n, extra;
   TYPE *C, *A, *B, *a, *b, *c;
   double t0, t1, mf;
   const TYPE alpha=1.0;
   TYPE beta=1.0;
   void *vp=NULL;

   if (NOMOVE)
   {
      ldc = (LDC) ? LDC : mb;
      setsz = (ldc * nb + kb*(mb+nb));
      vp = malloc(ATL_Cachelen + ATL_MulBySize(setsz));
      ATL_assert(vp);
      A =  ATL_AlignPtr(vp);
      B = A + mb*kb;
      C = B + kb*nb;
      for (i=0; i < setsz; i++) A[i] = dumb_rand();
      incA = incB = incC = 0;
   }
   else
   {
      if (movA && movB && movC)         /* no reuse at all */
      {
         setsz = ATL_MulBySize(mb*nb+kb*(mb+nb));
         nset = (FLSIZE+setsz-1)/setsz;
         FLSIZE = nset*setsz;
         setsz = mb*nb+kb*(mb+nb);
         vp = malloc(ATL_Cachelen + ATL_MulBySize(setsz));
         ATL_assert(vp);
         A = ATL_AlignPtr(vp);
         B = A + kb*mb*nset;
         C = B + kb*nb*nset;
         ldc = (LDC) ? mb*nset : mb;
         for (n=setsz*nset,i=0; i < n; i++) A[i] = dumb_rand();
         incA = mb*kb;
         incB = kb*nb;
         incC = mb*nb;
      }
      else if (movA && movB && !movC)   /* square-case ATLAS behavior */
      {
         setsz = kb*(mb+nb);
         ldc = (LDC) ? LDC : mb;
         ATL_assert(ldc >= mb);
         extra = ldc*nb;
         incA = mb*kb;
         incB = kb*nb;
         incC = 0;
      }
      else if (!movB && movA && movC)   /* rank-K behavior */
      {
         setsz = mb*(kb+nb);
         extra = kb*nb;
         incA = mb*kb;
         incB = 0;
         incC = mb*nb;
      }
      else
      {
         fprintf(stderr, "%s,%d: What case are you wanting?\n",
                 __FILE__, __LINE__);
         exit(-1);
      }
      if (!vp)
      {
         i = ATL_MulBySize(setsz);
         nset = (FLSIZE+i-1)/i;
         FLSIZE = nset * i;
         vp = malloc(ATL_Cachelen + ATL_MulBySize(FLSIZE+extra));
         ATL_assert(vp);
         A = ATL_AlignPtr(vp);
         if (movC)
         {
            C = A + mb*kb*nset;
            ldc = (LDC) ? mb*nset : mb;
            B = C + mb*nb*nset;
         }
         else
         {
            B = A + mb*kb*nset;
            C = B + kb*nb*nset;
         }
         for (n=setsz*nset+extra,i=0; i < n; i++) A[i] = dumb_rand();
      }
   }
   a = A; b = B; c = C;
   t0 = ATL_walltime();
   for (j=0,i=reps; i; i--)
   {
      #ifdef ATL_NEWTIME
         KMM(mblks, nblks, kb, a, b, c, movA ? a+incA : a,
             movB ? b+incB : b, movC ? c+incC : c);
      #else
         KMM(mb, nb, kb, alpha, a, kb, b, kb, beta, c, ldc);
      #endif
      if (++j != nset)
      {
         a += incA;
         b += incB;
         c += incC;
      }
      else
      {
         #ifndef ATL_NEWTIME
            beta = (beta != 0.0) ? -beta : 0.0;
         #endif
         j = 0;
         a = A; b = B; c = C;
      }
   }
   t1 = ATL_walltime() - t0;
   mf = (2.0*reps*mb*nb*kb) / (t1*1000000.0);
   free(vp);
   return(mf);
}

#define PRINT_COREID
#ifdef PRINT_COREID
   #include <utmpx.h>
#endif
//#define PRINT_NUMAIDS
#ifdef PRINT_NUMAIDS
   #define _GNU_SOURCE 1
   #include <unistd.h>
   #include <sys/syscall.h>
#endif
void *TimeOnCore(void *vp)
{
   struct kmm_struct *kp = vp;
   #ifdef PRINT_COREID
      printf("core=%d\n", sched_getcpu());
   #endif
#ifdef PRINT_NUMAIDS
    unsigned cpu, node;
    syscall(SYS_getcpu, &cpu, &node, NULL);
    printf("cpu=%u, node=%u\n", cpu, node);
#endif

   #ifdef ATL_NEWTIME
      kp->mf = GetKmmMflop(kp->mb, kp->nb, kp->kb, kp->mu, kp->nu, kp->ku,
                           kp->movA, kp->movB, kp->movC,
                           kp->FLSIZE, kp->reps, 0);
   #else
      kp->mf = GetKmmMflop(kp->mb, kp->nb, kp->kb, kp->movA, kp->movB, kp->movC,
                           kp->FLSIZE, kp->reps, kp->LDC);
   #endif
   return(NULL);
}


double *TimeOnCores(struct kmm_struct *kb)
{
   struct kmm_struct *kp;
   pthread_t *threads;
   pthread_attr_t *attr;
   cpu_set_t cpuset;
   double *mflops;
   int i, p;

   p = kb->p;
   kp = malloc(sizeof(struct kmm_struct)*p);
   threads = malloc(sizeof(pthread_t)*p);
   attr = malloc(sizeof(pthread_attr_t)*p);
   mflops = malloc(sizeof(double)*p);
   ATL_assert(kp && threads && attr && mflops);
   for (i=0; i < p; i++)
   {
      memcpy(kp+i, kb, sizeof(struct kmm_struct));
      kp[i].iam = i;
      CPU_ZERO(&cpuset);
      CPU_SET(kp->pids[i], &cpuset);
      assert(!pthread_attr_setaffinity_np(attr+i, sizeof(cpuset), &cpuset));
      pthread_create(threads+i, attr+i, TimeOnCore, kp+i);
   }
   for (i=0; i < p; i++)
   {
      pthread_join(threads[i], NULL);
      mflops[i] = kp[i].mf;
   }
   free(kp->pids);
   free(kp);
   free(threads);
   free(attr);
   return(mflops);
}

void GetStat(int n, double *d, double *min, double *max, double *avg)
{
   int i;
   double dmin, dmax, dsum;

   dmin = dmax = dsum = d[0];
   for (i=1; i < n; i++)
   {
      dmax = (dmax >= d[i]) ? dmax : d[i];
      dmin = (dmin <= d[i]) ? dmin : d[i];
      dsum += d[i];
   }
   *min = dmin;
   *max = dmax;
   *avg = dsum / (double)n;
}

void PrintUsage(char *name, int iarg, char *arg)
{
   fprintf(stderr, "\nERROR around arg %d (%s).\n", iarg, arg ? arg:"unknown");
   fprintf(stderr, "USAGE: %s [flags], where flags are:\n", name);
   fprintf(stderr, "   -p <#> : use # threads (with affinity)\n");
   fprintf(stderr, "   -P <#> id1 ... id# : spawn # threads to given IDs\n");
   fprintf(stderr, "   -B <#> : mb = nb = kb = #\n");
   fprintf(stderr, "   -m <#> : mb = #\n");
   fprintf(stderr, "   -n <#> : nb = #\n");
   fprintf(stderr, "   -k <#> : kb = #\n");
   #ifdef ATL_NEWTIME
   fprintf(stderr, "   -u[mnk] <#> : M/N/K loop unrolling is #\n");
   #endif
   fprintf(stderr, "   -r <#> : set the # of times to call KMM\n");
   fprintf(stderr, "   -F <kb> : set flush size in kilobytes\n");
   fprintf(stderr, "   -C <#> : set ldc; 0 means mb\n");
   fprintf(stderr, "   -M[a,b,c] <#> : mov[A,B,C] = #\n");
   exit(iarg ? iarg : -1);
}

struct kmm_struct *GetFlags(int nargs, char **args)
{
   struct kmm_struct *kp;
   int i, j;

   kp = malloc(sizeof(struct kmm_struct));
   ATL_assert(kp);
   kp->pids = NULL;
   kp->p = 1;
   kp->mb = kp->nb = kp->kb = 40;
   #ifdef ATL_NEWTIME
      kp->mu = kp->nu = 4;
      kp->ku = 1;
   #endif
   kp->movA = kp->movB = kp->movC = 0;
   kp->FLSIZE = L2SIZE;
   kp->reps = 200;
   kp->LDC = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'F':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->FLSIZE = atoi(args[i]) * 1024;
         break;
      case 'C':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->LDC = atoi(args[i]);
         break;
      #ifdef ATL_NEWTIME
      case 'u':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         j = atoi(args[i]);
         if (args[i-1][2] == 'k')
            kp->ku = j;
         else if (args[i-1][2] == 'n')
            kp->nu = j;
         else
            kp->mu = j;
         break;
      #endif
      case 'r':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->reps = atoi(args[i]);
         break;
      case 'P':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->p = atoi(args[i]);
         kp->pids = malloc(sizeof(int)*kp->p);
         assert(kp->pids);
         for (j=0; j < kp->p; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of arguments");
            kp->pids[j] = atoi(args[i]);
         }
         break;
      case 'p':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->p = atoi(args[i]);
         break;
      case 'm':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->mb = atoi(args[i]);
         break;
      case 'n':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->nb = atoi(args[i]);
         break;
      case 'k':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->kb = atoi(args[i]);
         break;
      case 'B':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->mb = kp->nb = kp->kb = atoi(args[i]);
         break;
      case 'M':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         switch(args[i-1][2])
         {
         case 'c':
         case 'C':
            kp->movC = atoi(args[i]);
            break;
         case 'b':
         case 'B':
            kp->movB = atoi(args[i]);
            break;
         case 'a':
         case 'A':
            kp->movA = atoi(args[i]);
            break;
         default:
            PrintUsage(args[0], i-1, "unknown mov matrix");
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!kp->pids)
   {
      kp->pids = malloc(sizeof(int)*kp->p);
      assert(kp->pids);
      #ifdef ATL_ARCH_XeonPHI
      {
         int n4 = ((kp->p)>>2)<<2, nr = kp->p - n4;
         for (j=0; j < n4; j += 2)
         {
            kp->pids[j] = 2*j;
            kp->pids[j+1] = 2*j+1;
         }
         switch(nr)
         {
         case 3:
            kp->pids[j+2] = 2*(j+2);
         case 2:
            kp->pids[j+1] = 2*j+1;
         case 1:
            kp->pids[j] = 2*j;
            break;
         case 0:;
         }
      }
      #else
         for (j=0; j < kp->p; j++)
             kp->pids[j] = j;
      #endif
   }
   return(kp);
}

int main(int nargs, char **args)
{
   struct kmm_struct *kp;
   int i, p;
   double *dp;
   double min, max, avg;
   FILE *fpout = stdout;

   kp = GetFlags(nargs, args);
   p = kp->p;
   dp = TimeOnCores(kp);
   free(kp);
   GetStat(p, dp, &min, &max, &avg);
   fprintf(fpout, "ALL CORES: min=%le, max=%le, avg=%le\n", min, max, avg);
   fprintf(fpout, "PER-CORE: %le", dp[0]);
   for (i=1; i < p; i++)
      fprintf(fpout, ", %le", dp[i]);
   fprintf(fpout, "\n\n%.2f\n", avg);
   free(dp);
   exit(0);
}
