#define ATL_TP_DECL 1
#include "atlas_taffinity.h"
#include "atlas_threads.h"
#include "atlas_misc.h"
#include "atlas_bitvec.h"

/*
 * Started wt pthread_create on affID==0, this guy will create the
 * rest of the work queue
 */
void *ATL_threadpool_launch(void *vp)
{
   ATL_thread_t *tp = vp;
   ATL_tpool_t *pp = ATL_TP_PTR;
   const int P = pp->nthr;
   int i;
/*
 * Set my affinity if I haven't already
 */
   #ifdef ATL_PAFF_SELF
      if (tp->affID < 0)
      {
         tp->affID = 1 - tp->affID;
         ATL_setmyaffinity(tp->affID);
      }
   #endif
   #ifndef ATL_TP_FULLPOLL
   if (!tp->rank)
      for (i=1; i < ATL_NTHREADS; i++)
         ATL_thread_start(tp+i, i, 0, ATL_threadpool_launch, tp+i);
   #endif

   ATL_threadpool(tp);
}
/*
 * This function is called only by one thread, and is unsafe if called
 * by more than one.  It is protected ty ATL_IsFirstThreadedCall, which
 * should return true for only one thread
 */
void ATL_InitThreadPoolStartup(int P, void *pd, void *extra)
{
   ATL_tpool_t *pp;
   int i;

   ATL_assert(!ATL_TP_PTR);


   ATL_TP_PTR = pp = ATL_NewThreadPool(ATL_NTHREADS, 1, pd);
   #if !defined(ATL_TP_FULLPOLL) && !defined(ATL_OMP_THREADS)
      ATL_mutex_lock(pp->tpmut);
   #endif
   pp->nworkers = P;
   pp->PD = pd;
   pp->extra = extra;
   #if defined(ATL_PHI_SORTED) || defined(ATL_PHI_SORTED2)
      ATL_TPF_SET_MICSORTED(pp);
   #endif

   #ifndef ATL_OMP_THREADS
      #ifdef ATL_TP_FULLPOLL
      for (i=1; i < ATL_NTHREADS; i++)
      {
         ATL_thread_start(pp->threads+i, i, 0, ATL_threadpool_launch,
                          pp->threads+i);
      }
      #else
         ATL_thread_start(pp->threads, i, 0, ATL_threadpool_launch,pp->threads);
      #endif
   #endif
/*
 * We leave function still holding the lock, so that master can go to sleep
 * awaiting job completion before spawned threads are able to complete
 * and thus mess up timing; ATL_TP_FULLPOLL never gets the lock, and will
 * have to call dojob() manually after allowing threads to proceed
 */
}

#if defined(ATL_TP_FULLPOLL) || defined(ATL_OMP_THREADS)
/*
 * This version uses cache to communicate with polling cores.
 */
void ATL_goParallel
(
   const unsigned int P0,  /* # of worker threads to use on job */
   void *DoWork,           /* ptr to function doing the job */
   void *DoComb,           /* ptr to func doing combine; NULL: don't combine */
   void *PD,               /* ptr to problem definition (work queue, etc) */
   void *extra             /* extra ptr that is passed to DoWork/DoComb */
)
{
   ATL_tpool_t *pp;
   int i, n, *ip;
   int JUST_STARTED=0;
   unsigned int P;
/*
 * Don't wakeup/create threads if being called for serial execution: just
 * do work on master process to save overhead & ease debugging
 */
   if (P0 < 2)

   {
      if (!ATL_TP_PTR)
         ATL_TP1_PTR = ATL_NewThreadPool(1, 0, NULL);
      pp = ATL_TP1_PTR;
      pp->nworkers = 1;
      pp->PD = PD;
      pp->extra = extra;
      pp->jobf = DoWork;
      pp->combf = DoComb;
      pp->jobf(pp, 0, 0);
      return;
   }
/*
 * If thread pool not currently active, must spawn it.
 */
   if (!ATL_TP_PTR)
   {
      #ifdef ATL_OMP_THREADS
      #pragma omp single
      #else
      if (ATL_IsFirstThreadedCall())
      #endif
      {
         JUST_STARTED = 1;
         ATL_setmyaffinity(0);
         ATL_InitThreadPoolStartup(Mmin(P0,ATL_NTHREADS), PD, extra);
      }
   }
   pp = ATL_TP_PTR;
   P = Mmin(P0, pp->nthr);
   pp->PD = PD;
   pp->extra = extra;
   pp->jobf = DoWork;
   pp->combf = DoComb;
   pp->nworkers = P;
/*
 * If I'm doing a combine, setup combine bitvectors based on nworkers
 */
   if (DoComb)
   {
      ATL_UnsetAllBitsBV(pp->combReadyBV);
      if (P == pp->nthr)
         ATL_UnsetAllBitsBV(pp->combDoneBV);
      else if (P >= pp->nthr - P) /* participating threads outnumber idle */
      {
         int i;
         ATL_UnsetAllBitsBV(pp->combDoneBV);
         for (i=P; i < pp->nthr; i++)
            ATL_SetBitBV(pp->combDoneBV, i);
      }
      else /* idle threads outnumber participating */
      {
         int i;
         ATL_SetAllBitsBV(pp->combDoneBV);
         for (i=0; i < P; i++)
            ATL_UnsetBitBV(pp->combDoneBV, i);
      }
   }
#ifdef ATL_OMP_THREADS
   omp_set_num_threads(P);
   #pragma omp parallel
   {
      int r;
      ATL_assert(omp_get_num_threads() == P);
      r = omp_get_thread_num();
      pp->threads[r].rank = r;
      pp->threads[r].P = P;
      #ifdef ATL_PAFF_SELF
         if (JUST_STARTED)
            ATL_setmyaffinity(r);
      #endif
      ATL_tpool_dojob(pp, r, 0);
   }
#else
/*
 * Tell workers to start working
 */
   for (i=1; i < P; i++)
   {
/*
 *    Wait for thread to get into busy-loop
 */
      while (pp->chkin[i] != 3)
         ATL_thread_yield();
/*
 *    Tell ready thread to start working
 */
      pp->chkin[i] = 1;
   }
   ATL_tpool_dojob(pp, 0, 0);
/*
 * Wait for all threads to complete
 */
   for (i=1; i < P; i++)
   {
      while (pp->chkin[i] != 3)
         ATL_thread_yield();
   }
#endif
}
#else
/*
 * Main way parallel jobs are spawned to workpool
 */
void ATL_goParallel
(
   const unsigned int P0,  /* # of worker threads to use on job */
   void *DoWork,           /* ptr to function doing the job */
   void *DoComb,           /* ptr to func doing combine; NULL: don't combine */
   void *PD,               /* ptr to problem definition (work queue, etc) */
   void *extra             /* extra ptr that is passed to DoWork/DoComb */
)
{
   ATL_tpool_t *pp;
   int n, *ip;
   int JUST_STARTED=0;
   unsigned int P;
/*
 * Don't wakeup/create threads if being called for serial execution: just
 * do work on master process to save overhead & ease debugging
 */
   if (P0 < 2)

   {
      if (!ATL_TP_PTR)
         ATL_TP1_PTR = ATL_NewThreadPool(1, 0, NULL);
      pp = ATL_TP1_PTR;
      pp->nworkers = 1;
      pp->PD = PD;
      pp->extra = extra;
      pp->jobf = DoWork;
      pp->combf = DoComb;
      pp->jobf(pp, 0, 0);
      return;
   }
/*
 * If thread pool not currently active, must spawn it.
 */
   if (!ATL_TP_PTR)
   {
      #ifdef ATL_OMP_THREADS
      #pragma omp single
      #else
      if (ATL_IsFirstThreadedCall())
      #endif
      {
         JUST_STARTED = 1;
         ATL_InitThreadPoolStartup(Mmin(P0,ATL_NTHREADS), PD, extra);
      }
   }
   pp = ATL_TP_PTR;
   P = Mmin(P0, pp->nthr);

   if (!JUST_STARTED)
   {
      ATL_mutex_lock(pp->tpmut);
      pp->jobID++;
      pp->PD = PD;
      pp->extra = extra;
      pp->jobf = DoWork;
      pp->combf = DoComb;
      pp->nworkers = P;
      pp->nwdone = pp->wcnt = 0;
      pp->WORKDONE = pp->NOWORK = 0;
      #ifdef ATL_TP_FORCEBCAST
         #if (ATL_TP_FORCEBCAST)
            ATL_TPF_UNSET_ZEROWAKES(pp);
         #else
            ATL_TPF_SET_ZEROWAKES(pp);
         #endif
      #else
         if (P == pp->nthr)
            ATL_TPF_UNSET_ZEROWAKES(pp);
         else
            ATL_TPF_SET_ZEROWAKES(pp);
         #endif
/*
 *    Wake threads up to do the new job
 */
      if (ATL_TPF_ZEROWAKES(pp))
         ATL_cond_signal(pp->wcond);
      else
      #ifdef ATL_PHI_SORTED
      {
         ATL_cond_bcast(pp->wcond);
         if (P > ATL_NTHREADS/4)
            ATL_cond_bcast(pp->wcond2);
         if (P > ATL_NTHREADS/2)
            ATL_cond_bcast(pp->wcond3);
         if (P > (3*ATL_NTHREADS)/4)
            ATL_cond_bcast(pp->wcond4);

      }
      #else
         ATL_cond_bcast(pp->wcond);
      #endif
   }
   else /* thread pool just started up, I have tpmut */
   {
      pp->jobf = DoWork;
      pp->combf = DoComb;
   }
/*
 * If I'm doing a combine, setup combine bitvectors based on nworkers
 */
   if (DoComb)
   {
      ATL_UnsetAllBitsBV(pp->combReadyBV);
      if (P == pp->nthr)
         ATL_UnsetAllBitsBV(pp->combDoneBV);
      else if (P >= pp->nthr - P) /* participating threads outnumber idle */
      {
         int i;
         ATL_UnsetAllBitsBV(pp->combDoneBV);
         for (i=P; i < pp->nthr; i++)
            ATL_SetBitBV(pp->combDoneBV, i);
      }
      else /* idle threads outnumber participating */
      {
         int i;
         ATL_SetAllBitsBV(pp->combDoneBV);
         for (i=0; i < P; i++)
            ATL_UnsetBitBV(pp->combDoneBV, i);
      }
   }
/*
 * Threads have started, so this routine awaits completion of task
 */
   do
   {
      #ifdef DEBUG
         fprintf(stderr, "master sleeps WD=%d, P=%d, comb=%p\n",
                 pp->WORKDONE, pp->nworkers, pp->combf);
      #endif
      ATL_cond_wait(pp->mcond, pp->tpmut);
      #ifdef DEBUG
         fprintf(stderr, "master awake  WD=%d\n", pp->WORKDONE);
      #endif
   }
   while (!pp->WORKDONE);
   #ifdef DEBUG
      fprintf(stderr, "MASTER AWAKE FINAL WD=%d\n", pp->WORKDONE);
   #endif
   pp->nworkers = 0;     /* make sure spurious wakeup will go back to sleep */
/*
 * If I'm using a thread pool where all threads are known to be asleep before
 * returning, then on the first launch we must make sure they have all checked
 * in before returning (in subsequent uses of pool, master will not be
 * awakened until all workers have gone to sleep)
 */
   #ifndef ATL_POLLTPOOL
      if (JUST_STARTED && pp->nsleep < pp->nthr)
      {
         do
         {
            ATL_mutex_unlock(pp->tpmut);
            while(pp->nsleep < pp->nthr)
               ATL_thread_yield();
            ATL_mutex_lock(pp->tpmut);
         }
         while(pp->nsleep < pp->nthr);
      }
   #endif
   ATL_mutex_unlock(pp->tpmut);
}
#endif

ATL_tpool_t *ATL_NewThreadPool
(
   const int P,  /* length of thread array to allocate, 0 for don't */
   int ICOM,     /* 0: don't allocate icomm array */
   void *vp      /* what to initialize thread arrays vp to */
)
{
   ATL_tpool_t *pp;
   pp = calloc(1, sizeof(ATL_tpool_t));
   ATL_assert(pp);
   if (P)
   {
      int i;
      pp->nthr = P;
      pp->threads = malloc(P*sizeof(ATL_thread_t));
      ATL_assert(pp->threads);
      if (!vp)
         vp = pp;
      for (i=0; i < P; i++)
      {
         pp->threads[i].rank = i;
         pp->threads[i].P = P;
         pp->threads[i].affID = -i-1;
         pp->threads[i].vp = vp;
      }
   }
   #ifdef ATL_TP_FULLPOLL
      pp->chkin = calloc(P, sizeof(short));
   #elif !defined(ATL_OMP_THREADS) && !defined(ATL_WINTHREADS)
      pp->jobID = 0;
      pp->mcond = ATL_cond_init();
      pp->wcond = ATL_cond_init();
   #endif
   if (ICOM)
   {
      pp->icomm = calloc(P, sizeof(int));
      ATL_assert(pp->threads);
   }
   pp->combmut = ATL_mutex_init();
   pp->tpmut = ATL_mutex_init();
   #ifdef ATL_PHI_SORTED
      pp->wcond2 = ATL_cond_init();
      pp->wcond3 = ATL_cond_init();
      pp->wcond4 = ATL_cond_init();
   #elif defined(ATL_PHI_SORTED2)
   {
      int i, n;
      pp->ncntxts = 0;
      n = ATL_NTHREADS >> 2;
      pp->wconds = malloc(sizeof(void*)*n);
      ATL_assert(pp->wconds);
      for (i=0; i < n; i++)
         pp->wconds[i] = ATL_cond_init();
   }
   #endif
   pp->combReadyBV = ATL_NewBV(P);
   pp->combDoneBV = ATL_NewBV(P);
   return(pp);
}

void ATL_FreeThreadPool(ATL_tpool_t *pp)
{
   if (!pp)
      return;
   if (pp->icomm)
      free((void*)pp->icomm);
   if (pp->threads)
      free(pp->threads);
   if (pp->combmut)
      ATL_mutex_free(pp->combmut);
   if (pp->tpmut)
      ATL_mutex_free(pp->tpmut);
   #ifdef ATL_TP_FULLPOLL
      free(pp->chkin);
   #elif !defined(ATL_OMP_THREADS) && !defined(ATL_WINTHREADS)
   if (pp->mcond)
      ATL_cond_free(pp->mcond);
   if (pp->wcond)
      ATL_cond_free(pp->wcond);
   #endif
   #if defined(ATL_PHI_SORTED)
      if (pp->mcond2)
         ATL_cond_free(pp->mcond2);
      if (pp->mcond3)
         ATL_cond_free(pp->mcond3);
   #elif defined(ATL_PHI_SORTED2)
   {
      int i, n = pp->nthr>>2;
      for (i=0; i < n; i++)
         ATL_cond_free(pp->mconds[i]);
      free(pp->mconds);
   }
   #endif
   if (pp->combReadyBV)
      ATL_FreeBV(pp->combReadyBV);
   if (pp->combDoneBV)
      ATL_FreeBV(pp->combDoneBV);
   free(pp);
}
