#define _GNU_SOURCE 1 /* what manpage says you need to get CPU_SET */
#define __USE_GNU   1 /* what you actually have to set on linuxes I've seen */
#include <sched.h>    /* must include this before pthreads */
#include <pthread.h>
#include "atlas_misc.h"
#include <assert.h>
#include <string.h>
#define dumb_prand(is_) ( 0.5 - ((double)rand_r(is_))/((double)RAND_MAX) )


#define CINT const int

#ifdef TCPLX
   size_t rszA, rszB, rszC;
   void CAMM_b0(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
                const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
                const TYPE *pCn);
   void CAMM_b1(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
                const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
                const TYPE *pCn);
   void CAMM_bn(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
                 const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
                 const TYPE *pCn);
   #ifdef BETA0
      void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *Ai,
               const TYPE *Bi, TYPE *Ci, const TYPE *pAn, const TYPE *pBn,
               const TYPE *pCn)
      {
         extern size_t rszA, rszB, rszC;
         const TYPE *Ar=Ai+rszA, *Br=Bi+rszB;
         TYPE *Cr = Ci + rszC;
         CAMM_b0(mblks, nblks, K, Ai, Bi, Cr, Ai, Br, Ci);
         CAMM_b0(mblks, nblks, K, Ai, Br, Ci, Ar, Br, Cr);
         CAMM_bn(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
         CAMM_b1(mblks, nblks, K, Ar, Bi, Ci, pAn, pBn, pCn);
      }
   #elif defined(BETA1)
      void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *Ai,
               const TYPE *Bi, TYPE *Ci, const TYPE *pAn, const TYPE *pBn,
               const TYPE *pCn)
      {
         extern size_t rszA, rszB, rszC;
         const TYPE *Ar=Ai+rszA, *Br=Bi+rszB;
         TYPE *Cr = Ci + rszC;
         CAMM_bn(mblks, nblks, K, Ai, Bi, Cr, Ai, Br, Ci);
         CAMM_b1(mblks, nblks, K, Ai, Br, Ci, Ar, Br, Cr);
         CAMM_bn(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
         CAMM_b1(mblks, nblks, K, Ar, Bi, Ci, pAn, pBn, pCn);
      }
   #else
      void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *Ai,
               const TYPE *Bi, TYPE *Ci, const TYPE *pAn, const TYPE *pBn,
               const TYPE *pCn)
      {
         extern size_t rszA, rszB, rszC;
         const TYPE *Ar=Ai+rszA, *Br=Bi+rszB;
         TYPE *Cr = Ci + rszC;
         CAMM_b1(mblks, nblks, K, Ai, Bi, Cr, Ai, Br, Ci);
         CAMM_bn(mblks, nblks, K, Ai, Br, Ci, Ar, Br, Cr);
         CAMM_b1(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
         CAMM_bn(mblks, nblks, K, Ar, Bi, Ci, pAn, pBn, pCn);
      }
   #endif
#else
   #ifndef KMM
      #define KMM ATL_USERMM
   #endif
   void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
            const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
            const TYPE *pCn);
#endif

struct kmm_struct{
   int mb, nb, kb;                      /* C: mbxnb, At: kbxmb, B: kbXnb */
   int mu, nu, ku;                      /* needed to compute mblks/nblks */
   int movA, movB, movC;                /* which mat move in flush array? */
   int FLSIZE;                          /* min area to move in in bytes */
   int reps;                            /* # calls to kmm in one timing */
   int LDC;                             /* what should ldc be set to? */
   int iam;                             /* thread rank */
   int p;                               /* total number of threads */
   int *pids;                           /* IDs of processors for affinity */
   double mf;                           /* mflop returned by timing */
   volatile unsigned char *chkin;       /* P-len array to signal start/done */
};

double GetKmmMflop
(
   CINT mb, CINT nb, CINT kb,           /* C: mbxnb, At: kbxmb, B: kbXnb */
   CINT mu, CINT nu, CINT ku,
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
   CINT mblks = mb/mu, nblks = nb/nu;
   const int NOMOVE = !(movA|movB|movC);
   size_t ldc, setsz, nset, i, j, incA, incB, incC, n, extra, EXTRA=0;
   TYPE *C, *A, *B, *a, *b, *c;
   double t0, t1, mf;
   const TYPE alpha=1.0;
   TYPE beta=1.0;
   void *vp=NULL;
   unsigned int seed = mb*kb + (nb<<14);
   EXTRA = 4*mu*nu;


   if (NOMOVE)
   {
      setsz = (mb * nb + kb*(mb+nb));
      vp = malloc(2*ATL_Cachelen + ATL_MulBySize(setsz+EXTRA));
      ATL_assert(vp);
      A = ATL_AlignPtr(vp);
      B = A + (mb*kb SHIFT);
      C = B + (kb*nb SHIFT);
      C = ATL_AlignPtr(C);
      #ifdef TCPLX
         setsz += setsz;
      #endif
      for (i=0; i < setsz; i++) A[i] = dumb_prand(&seed);
      incA = incB = incC = 0;
   }
   else
   {
      if (movA && movB && movC)         /* no reuse at all */
      {
         setsz = ATL_MulBySize(mb*nb+kb*(mb+nb));
         nset = (FLSIZE+setsz-1)/setsz;
         FLSIZE = nset*(setsz+ATL_Cachelen);
         setsz = mb*nb + kb*(mb+nb) + ATL_DivBySize(ATL_Cachelen);
         vp = malloc(ATL_Cachelen + ATL_MulBySize(nset*setsz+EXTRA));
         ATL_assert(vp);
         setsz = setsz SHIFT;
         A = ATL_AlignPtr(vp);
         B = A + (kb*mb*nset SHIFT);
         C = B + (kb*nb*nset SHIFT);
         C = ATL_AlignPtr(B);
         for (n=setsz*nset,i=0; i < n; i++) A[i] = dumb_prand(&seed);
         incA = mb*kb SHIFT;
         incB = kb*nb SHIFT;
         incC = mb*nb SHIFT;
      }
      else if (movA && movB && !movC)   /* square-case ATLAS behavior */
      {
         setsz = kb*(mb+nb);
         extra = mb*nb;
         incA = mb*kb SHIFT;
         incB = kb*nb SHIFT;
         incC = 0;
      }
      else if (!movB && movA && movC)   /* rank-K behavior */
      {
         setsz = mb*(kb+nb);
         extra = kb*nb;
         incA = mb*kb SHIFT;
         incB = 0;
         incC = mb*nb SHIFT;
      }
      else if (movA && !movB && !movC) /* used to det cost of B/C access */
      {
         setsz = mb*kb;
         incA = setsz SHIFT;
         incB = incC = 0;
         extra = mb*nb + kb*nb;
      }
      else if (!movA && movB && !movC) /* used to det cost of ld A/C access */
      {
         incB = nb*kb SHIFT;
         incA = incC = 0;
         setsz = incB;
         extra = mb*(nb+kb);
      }
      else if (!movA && !movB && movC) /* used to det cost of ld A/B access */
      {
         incC = mb*nb SHIFT;
         incA = incB = 0;
         setsz = incC;
         extra = kb*(mb+nb);
      }
      else
      {
         fprintf(stderr, "%s,%d: What case are you wanting?\n",
                 __FILE__, __LINE__);
         exit(-1);
      }
      #ifdef TCPLX
         rszA = mb*kb;
         rszB = kb*nb;
         rszC = mb*nb;
      #endif
      if (!vp)
      {
         TYPE *dp, *ep;
         i = ATL_MulBySize(setsz);
         nset = (FLSIZE+i-1)/i;
         FLSIZE = nset * i;
         vp = malloc(6*ATL_Cachelen + ATL_MulBySize(extra+EXTRA) + FLSIZE);
         ATL_assert(vp);
         ep = ATL_AlignPtr(vp);
         dp = (TYPE*)(((char*)ep) + 2*ATL_Cachelen + ATL_MulBySize(extra));
         if (movA)
         {
            A = ATL_AlignPtr(dp);
            dp = A + incA*nset;
         }
         else
         {
            A = ep;
            ep += (mb*kb SHIFT);
         }
         if (movB)
         {
            B = ATL_AlignPtr(dp);
            dp = B + incB*nset;
         }
         else
         {
            B = ATL_AlignPtr(ep);
            ep = B + (kb*nb SHIFT);
         }
         C = (movC) ? dp : ep;
         C = ATL_AlignPtr(C);
         dp = ATL_AlignPtr(vp);
         setsz = setsz SHIFT;
         n = setsz*nset + (extra SHIFT);
         for (i=0; i < n; i++)
            dp[i] = dumb_prand(&seed);
      }
   }
   a = A; b = B; c = C;
   t0 = ATL_walltime();
   for (j=0,i=reps; i; i--)
   {
      KMM(mblks, nblks, kb, a, b, c, movA ? a+incA : a,
          movB ? b+incB : b, movC ? c+incC : c);
      if (++j != nset)
      {
         a += incA;
         b += incB;
         c += incC;
      }
      else
      {
         j = 0;
         a = A; b = B; c = C;
      }
   }
   t1 = ATL_walltime() - t0;
   mf = (2.0*reps*mb*nb*kb) / (t1*1000000.0);
   #ifdef TCPLX
      mf *= 4.0;
   #endif
   free(vp);
   return(mf);
}

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
   const int P = kp->p;
   int i;
   volatile unsigned char *chkin = kp->chkin;
   #ifdef PRINT_COREID
      printf("core=%d\n", sched_getcpu());
   #endif
#ifdef PRINT_NUMAIDS
    unsigned cpu, node;
    syscall(SYS_getcpu, &cpu, &node, NULL);
    printf("cpu=%u, node=%u\n", cpu, node);
#endif
/*
 * First we barrier, so that all cores are active.  Otherwise, 1st core to
 * start may run with bus to himself, and not give us a measure of true
 * parallel performance.  Want full-on contention as in perfect parallel code.
 * We wait on chkin array to have all non-zero entries.  Even on weakly-ordered
 * caches this should work, though the delay may be long.
 */
   chkin[kp->iam] = 1;
   for (i=0; i < P; i++)
      while(!chkin[i]);
   kp->mf = GetKmmMflop(kp->mb, kp->nb, kp->kb, kp->mu, kp->nu, kp->ku,
                        kp->movA, kp->movB, kp->movC, kp->FLSIZE, kp->reps, 0);
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
   unsigned char *chkin;

   p = kb->p;
   kp = malloc(sizeof(struct kmm_struct)*p);
   threads = malloc(sizeof(pthread_t)*p);
   attr = malloc(sizeof(pthread_attr_t)*p);
   mflops = malloc(sizeof(double)*p);
   chkin = malloc(sizeof(char)*p);
   ATL_assert(kp && threads && attr && mflops && chkin);
   for (i=0; i < p; i++)   /* init chkin to 0 before starting any threads */
      chkin[i] = 0;        /* when all entries non-zero, all thrds started */
   for (i=0; i < p; i++)
   {
      memcpy(kp+i, kb, sizeof(struct kmm_struct));
      kp[i].chkin = (volatile char*)chkin;
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
   fprintf(stderr, "   -u[mnk] <#> : M/N/K loop unrolling is #\n");
   fprintf(stderr, "   -r <#> : set the # of times to call KMM\n");
   fprintf(stderr, "   -R <mf>: set # reps to force <mf> MFLOPs\n");
   fprintf(stderr, "   -F <kb> : set flush size in kilobytes\n");
   fprintf(stderr, "   -M[a,b,c] <#> : mov[A,B,C] = #\n");
   exit(iarg ? iarg : -1);
}

struct kmm_struct *GetFlags(int nargs, char **args, FILE **fpout)
{
   struct kmm_struct *kp;
   double mflops=750.0;
   int i, j;

   *fpout = NULL;
   kp = malloc(sizeof(struct kmm_struct));
   ATL_assert(kp);
   kp->pids = NULL;
   kp->p = 1;
   kp->mb = kp->nb = kp->kb = 40;
   kp->mu = kp->nu = 4;
   kp->ku = 1;
   kp->movA = kp->movB = kp->movC = 0;
   kp->FLSIZE = L2SIZE;
   kp->reps = 0;
   kp->LDC = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'f':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *fpout = fopen(args[i], "w");
         break;
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
      case 'R':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         mflops = atof(args[i]);
         break;
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
   if (!kp->reps)
   {
      kp->reps = (mflops*1000000.0/((2.0*kp->mb)*kp->nb*kp->kb));
      if (kp->reps < 1)
         kp->reps = 1;
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

   kp = GetFlags(nargs, args, &fpout);
   p = kp->p;
   dp = TimeOnCores(kp);
   free(kp);
   GetStat(p, dp, &min, &max, &avg);
   printf("PER-CORE: %le", dp[0]);
   for (i=1; i < p; i++)
      printf(", %le", dp[i]);
   printf("\nALL CORES: min=%.2f, max=%.2f, avg=%.2f\n", min, max, avg);
   if (fpout)
   {
      fprintf(fpout, "%d 1\n", p);
      for (i=0; i < p; i++)
         fprintf(fpout, "%e\n", dp[i]);
   }
   free(dp);
   if (fpout && fpout != stdout && fpout != stderr)
      fclose(fpout);
   exit(0);
}
