#include "atlas_threads.h"
#include "atlas_misc.h"
#include "atlas_bitvec.h"
/*
 * These are wrapperf for original codes, which expect to get a launchstruct
 * as the problem definition, and DoWork expects to get a pointer to my entry
 * in the thread array
 */
void ATL_oldjobwrap(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp = vpp;
   ATL_LAUNCHSTRUCT_t *lp = pp->PD;
   void *vp = pp->threads[rank].vp;

   pp->threads[rank].vp = lp;
   lp->DoWork(lp, pp->threads+rank);
   pp->threads[rank].vp = vp;
}

void ATL_oldcombwrap(void *vpp, int rank, int vrank, int vhisrank)
{
   ATL_tpool_t *pp = vpp;
   ATL_LAUNCHSTRUCT_t *lp = pp->PD;

   if (lp->DoComb)
      lp->DoComb(lp->opstruct, vrank, vhisrank);
}

/*
 * Fully dynamic combine with restriction that iam=0 gets final answer
 * NOTE: assumed that during operation, have change tp->P and tp->rank
 *       to local values.  They will be restored to true rank & P here.
 */
static void ATL_dyncomb(ATL_tpool_t *pp, const int rank, const int vrank)
{
   const unsigned int P=pp->nworkers;
   ATL_thread_t *tp = pp->threads + rank;
/*
 * Let everyone know my data is ready for combining
 * Thread 0's work is not available because some combines expect
 * thread 0 to have the answer
 */
   if (vrank)
   {
      ATL_mutex_lock(pp->combmut);
      ATL_SetBitBV(pp->combReadyBV, vrank);
      ATL_mutex_unlock(pp->combmut);
   }
/*
 * Now participate in combine until everything is done, or someone
 * combines my stuff, and thus ends my participation
 */
   do
   {
      int d;
/*
 *    If my combDone bit is set, someone has combined my results: I'm finished
 */
      if (ATL_IsBitSetBV(pp->combDoneBV, vrank))
         return;
/*
 *    Seize mutex lock and do the work if still there & I'm still in game
 */
      ATL_mutex_lock(pp->combmut);
      if (ATL_IsBitSetBV(pp->combDoneBV, vrank)) /* I'm done */
      {
         ATL_mutex_unlock(pp->combmut);
         return;
      }
      d = ATL_FindFirstSetBitBV(pp->combReadyBV, 0);
      if (d == vrank)
         d = (vrank != P-1) ?
            ATL_FindFirstSetBitBV(pp->combReadyBV, vrank+1) : -1;
      if (d != -1)  /* work still needs to be done */
      {
         ATL_SetBitBV(pp->combDoneBV, d);        /* combined thr is done */
         ATL_UnsetBitBV(pp->combReadyBV, d);     /* I took, no longer avail */
         ATL_UnsetBitBV(pp->combReadyBV, vrank); /* I'm not avail, working */
         ATL_mutex_unlock(pp->combmut);
         pp->combf(pp, rank, vrank, d);
/*
 *       After finishing my combine, my buffer ready for other's use
 */
         if (vrank)  /* thr 0 must get answer, so never available */
         {
            ATL_mutex_lock(pp->combmut);
            ATL_SetBitBV(pp->combReadyBV, vrank);   /* now avail again */
         }
/*
 *       Check if I'm the last guy left, if so, combine is complete
 *       Only vrank=0 allowed to be last guy (0's buff has final answer)
 */
         else  /* I'm 0, so see if combine is done */
         {
            ATL_mutex_lock(pp->combmut);
            if (ATL_FindFirstUnsetBitBV(pp->combDoneBV, 1) == -1)
            {
               ATL_SetBitBV(pp->combDoneBV, vrank);    /* I'm done */
               ATL_UnsetBitBV(pp->combReadyBV, vrank); /* so not ready */
               ATL_mutex_unlock(pp->combmut);
               return;
            }
         }
      }
      else if (vrank)  /* No work to do, leave to reduce lock contention */
      {
         ATL_mutex_unlock(pp->combmut);
         return;
      }
      ATL_mutex_unlock(pp->combmut);
   }
   while(1);   /* end combine loop */
}

/*
 * This combine is for routines that require only leftward combination,
 * so that you can add results only from vranks greater (to the right)
 * than yourself, and then only if you have already added in all the
 * ranks in between yourself and the candidate.  This type of result
 * is naturally enforced by log2 reduction and linear sum-to-zero, which
 * some of the earlier tblas assume.
 */
static void ATL_leftcomb(ATL_tpool_t *pp, const int rank, const int vrank)
{
   const unsigned P = pp->nthr;
   const unsigned int nw = pp->nworkers;
   ATL_thread_t *tp = pp->threads + rank;
/*
 * Let everyone know my data is ready for combining
 * Thread 0's work is not available because some combines expect
 * thread 0 to have the answer
 */
   if (vrank)
   {
      ATL_mutex_lock(pp->combmut);
      ATL_SetBitBV(pp->combReadyBV, vrank);
      ATL_mutex_unlock(pp->combmut);
   }
/*
 * Now participate in combine until everything is done, or someone
 * combines my stuff, and thus ends my participation
 */
   do
   {
      int i, d;
/*
 *    If my combDone bit is set, somone has combined my results, so I'm finished
 */
      if (ATL_IsBitSetBV(pp->combDoneBV, vrank))
         return;
/*
 *    If everyone to right of me is finished, there's nothing left for me to do
 */
      for (i=vrank+1; i < nw && ATL_IsBitSetBV(pp->combDoneBV, i); i++);
      if (i == nw)
         return;
/*
 *    Sieze mutex lock and do the work if I'm still in the game
 */
      ATL_mutex_lock(pp->combmut);
      if (ATL_IsBitSetBV(pp->combDoneBV, vrank)) /* I'm done */
      {
         ATL_mutex_unlock(pp->combmut);
         return;
      }
/*
 *    Take buffers from as far right as possible, as they have least to do
 */
      for (i=nw-1; i > vrank; i--)
      {
/*
 *       If target is ready, there must be no unfinished nodes between us
 */
         if (ATL_IsBitSetBV(pp->combReadyBV, i)) /* got a candidate */
         {
            int j;
            for (j=vrank+1; j < i; j++)
               if (!ATL_IsBitSetBV(pp->combDoneBV, j))
                  break;
            if (j == i)  /* he's a legal partner */
               break;
         }
      }
      if (i != vrank)  /* work still needs to be done */
      {
         ATL_SetBitBV(pp->combDoneBV, i);        /* he's done, can quit */
         ATL_UnsetBitBV(pp->combReadyBV, i);     /* he's not avail, done  */
         ATL_UnsetBitBV(pp->combReadyBV, vrank); /* I'm not avail, working */
         ATL_mutex_unlock(pp->combmut);
         pp->combf(pp, rank, vrank, i);
/*
 *       After finishing my combine, my buffer ready for other's use
 */
         if (vrank)    /* rank=0 gets answer, so never available */
         {
            ATL_mutex_lock(pp->combmut);
            ATL_SetBitBV(pp->combReadyBV, vrank);   /* now avail again */
         }
/*
 *       Check if I'm the last guy left, if so, combine is complete
 *       Only vrank=0 allowed to be last guy (0's buff has final answer)
 */
         else  /* I'm 0, so see if combine is done */
         {
            ATL_mutex_lock(pp->combmut);
            if (ATL_FindFirstUnsetBitBV(pp->combDoneBV, 1) == -1)
            {
               tp->rank = rank;                       /* put my thr info back */
               tp->P = P;                             /* to global values */
               ATL_SetBitBV(pp->combDoneBV, vrank);   /* I'm done */
               ATL_UnsetBitBV(pp->combReadyBV, vrank);/* so not ready */
               ATL_mutex_unlock(pp->combmut);
               return;
            }
         }                             /* end else I'm 0 */
      }                                /* end work still to be done */
      ATL_mutex_unlock(pp->combmut);
   }
   while(1);   /* end combine loop */
}

#ifdef ATL_OMP_THREADS
/*
 * Do a task using thread pool pp; RETURNS: virtual rank (0...nworkers-1)
 */
int ATL_tpool_dojob
(
   ATL_tpool_t *pp,    /* thread pool to use */
   const int rank,     /* actual rank, my tp found at pp->threads+rank */
   const int CFWTHR    /* unused in OpenMP version */
)
{
/*
 * Do work based on PD
 */
   pp->jobf(pp, rank, rank);     /* do the required task */
/*
 * If the calling routine has requested it, perform dynamic combine
 */
   if (pp->combf)
   {
/*
 *    Do combine using virtial ranks stored in pp->icomm
 */
      if (ATL_TPF_DYNCOMB(pp))
         ATL_dyncomb(pp, rank, rank);
      else
         ATL_leftcomb(pp, rank, rank);
   }     /* end combine if */

   return(rank);
}

void *ATL_threadpool(void *vp)
{
   ATL_assert(0);
}
#elif !defined(ATL_TP_FULLPOLL)
/*
 * Do a task using thread pool pp; RETURNS: virtual rank (0...nworkers-1)
 */
int ATL_tpool_dojob
(
   ATL_tpool_t *pp,    /* thread pool to use */
   const int rank,     /* actual rank, my tp found at pp->threads+rank */
   const int CFWTHR    /* 1: worker thr called, 0: ATL_threadpool loop called */
)
{
   const int nworkers=pp->nworkers, nthr = pp->nthr;
   int vrank, iwrk;
   void *tpmut = pp->tpmut;
   ATL_thread_t *tp = pp->threads+rank;
/*
 * Enroll thread in working pool
 */
   if (CFWTHR)
      ATL_mutex_lock(tpmut);
/*
 * On PHI, we must use our actual ranks when we are exploiting contexts
 */
   #if defined(ATL_PHI_SORTED)
   if (nworkers == nthr || nworkers == (nthr>>1) ||
       nworkers == (nthr>>2)+(nthr>>1))
      vrank = rank;
   else
/*
 * For PHI_SORTED2, we always use all contexts, and only context 3 is awoken
 * by the usual process.  Context 3 then awakens contexts 0-2 on its node.
 */
   #elif defined(ATL_PHI_SORTED2)
      if (pp->ncntxts > 1)
      {
         volatile int *chkin = (volatile int*)
                      (((volatile char *)pd->chkin)+ATL_Cachelen*icore);
/*
 *       context 3 increments vrank count by ncontexts
 */
         if (rank & 3 == 3)
            chkin[3] = -pp->wcnt;
         else
            while (chkin[3] < ATL_NTHREADS);
         vrank = (rank&3) - chkin[3];
      }
      else
   #endif
      vrank = pp->wcnt;
   if (vrank >= nworkers || pp->wcnt >= nworkers) /* If I'm not a worker */
   {
      ATL_mutex_unlock(tpmut);                    /* release mutex and */
      return(vrank);                              /* get out of here now */
   }
   #if defined(ATL_PHI_SORTED2)
      if (pp->ncntxts > 1)
      {
         iwrk = pp->wcnt;
         pp->wcnt += pp->ncntxts;
      }
      else
   #endif
   iwrk = pp->wcnt++;        /* inc count of participating workers */
   ATL_mutex_unlock(tpmut);  /* let other threads enroll in job */
/*
 * If this was called from threadloop, I may need to wake other workers
 * Only if I'm the first worker to get here (!iwrk) & called from tl (!CFWTHR)
 */
   if (!(CFWTHR+iwrk) && ATL_TPF_ZEROWAKES(pp))
   {
   #ifdef ATL_PHI_SORTED
      if (ATL_TPF_MICSORTED(pp))
      {
         int k, nw=nworkers;
         const int p4 = (pp->nthr >> 2);

         if (nworkers < p4)
            for (k=1; k < p4; k++)
               ATL_cond_signal(pp->wcond);
         else
         {
            ATL_cond_bcast(pp->wcond);
            nw -= p4;
            if (nw < p4)
               while (nw--)
                  ATL_cond_signal(pp->wcond2);
            else
            {
               ATL_cond_bcast(pp->wcond2);
               nw -= p4;
               if (nw < p4)
                  while (nw--)
                     ATL_cond_signal(pp->wcond3);
               else
               {
                  ATL_cond_bcast(pp->wcond3);
                  nw -= p4;
                  if (nw <  p4)
                     while (nw--)
                        ATL_cond_signal(pp->wcond4);
                  else
                     ATL_cond_bcast(pp->wcond4);
               }
            }
         }
      }
      else
   #endif
      if (nthr == nworkers)
         ATL_cond_bcast(pp->wcond);
      else
      {
         int k;
         #ifdef ATL_PHI_SORTED2
         if (pp->ncntxts > 1)
         {
            for (k=nworkers-4; k > 0; k -= 4)
               ATL_cond_signal(pp->wcond);
         }
         else
         #endif
         for (k=1; k < nworkers; k++)
            ATL_cond_signal(pp->wcond);
      }
   }
/*
 * Wake up rest of threads on my core if contexts are called for
 */
   #ifdef ATL_PHI_SORTED2
      if (rank & 3 == 3 && pp->ncntxts > 1)
         ATL_cond_bcast(pp->wconds[icore]);
   #endif
/*
 * Now, do work based on PD; change thread info to use virtual rank info
 */
   tp->P = nworkers;             /* make jobf think only nworkers thrs */
   tp->rank = vrank;             /* and use virtual rank */
   pp->jobf(pp, rank, vrank);    /* do the required task */
/*
 * If the calling routine has requested it, perform dynamic combine
 */
   if (pp->combf)
   {
/*
 *    Do combine using virtial ranks stored in pp->icomm
 */
      if (ATL_TPF_DYNCOMB(pp))
         ATL_dyncomb(pp, rank, vrank);
      else
         ATL_leftcomb(pp, rank, vrank);
   }     /* end combine if */
/*
 * Return thread values to their global values
 */
   tp->rank = rank;
   tp->P = nthr;

   return(vrank);
}
#endif

#if defined(ATL_TP_FULLPOLL) && !defined(ATL_OMP_THREADS)
/*
 * Do a task using thread pool pp; RETURNS: virtual rank (0...ncores*4-1)
 */
int ATL_tpool_dojob
(
   ATL_tpool_t *pp,    /* thread pool to use */
   const int rank,     /* actual rank, my tp found at pp->threads+rank */
   const int CFWTHR    /* 1: worker thr called, 0: ATL_threadpool loop called */
)
{
   ATL_thread_t *tp = pp->threads+rank;
   const int nworkers=pp->nworkers, nthr = pp->nthr;
   #ifdef ATL_PHI_SORTED
      const int mycont = rank / ncores;
      int vrank, iwrk;
   #else
      #define vrank rank
   #endif

/*
 * Map ranks sorted by cores (0...ncores-1, ncores...ncores*2-1,...) to one
 * where all 4 cores are sequential
 */

   #if ATL_PHI_SORTED
      if (pp->ncntxts > 1)
      {
         vrank = (coreid<<2) + mycont;
         if (coreid > ncores || mycont > pp->ncntxts)
            return(vrank);
         tp->P = nworkers<<2;          /* make jobf think only nworkers thrs */
         tp->rank = vrank;             /* and use virtual rank */
      }
      else
      {
         if (rank >= nworkers)
            return(rank);
         vrank = rank;
         tp->P = nworkers;
      }
   #else
      if (rank >= nworkers)
         return;
      tp->P = nworkers;
   #endif
/*
 * Now, do work based on PD; change thread info to use virtual rank info
 */
   pp->jobf(pp, rank, vrank);    /* do the required task */
/*
 * If the calling routine has requested it, perform dynamic combine
 */
   if (pp->combf)
   {
/*
 *    Do combine using virtial ranks stored in pp->icomm
 */
      if (ATL_TPF_DYNCOMB(pp))
         ATL_dyncomb(pp, rank, vrank);
      else
         ATL_leftcomb(pp, rank, vrank);
   }     /* end combine if */
/*
 * Return thread values to their global values
 */
   #ifdef ATL_PHI_SORTED
      tp->rank = rank;
   #endif
   tp->P = nthr;

   return(vrank);
}
   #ifndef ATL_PHI_SORTED
      #undef vrank
   #endif
/*
 * This version of threadpool expects master process is bound to a non-used
 * core, and all workers area always polling on volatile variables for whether
 * to work.  Since cache coherence is orders of magnitude faster than cond
 * vars or mutexes, this should always win unless it causes contention with
 * non-ATLAS threads.  For instance, on a 56-core PHI, the last thread to wake
 * up using cond vars (mutexes) is 2835 (87) times slower than waking up
 * with cache coherence!
 */
void *ATL_threadpool(void *vp)
{
   ATL_thread_t *tp=vp;
   ATL_tpool_t *pp = ATL_TP_PTR;
   const unsigned int rank=tp->rank;
   const unsigned int nthr=pp->nthr;

   #ifdef DEBUG
      fprintf(stderr, "%d: WORKER START P=%d\n", rank, pp->nthr);
   #endif
/*
 * Infinite loop over tasks.  Loop has 2 phases:
 * (1) Poll for new work
 * (2) Do work if commanded, including optional combine
 */
   while (1)
   {
/*
 *    Wait for work to be posted.
 */
      pp->chkin[rank] = 3;  /* tell master in polling loop */
      do
      {
         ATL_thread_yield();
      }
      while (pp->chkin[rank] != 1);  /* await signal to do work */
/*
 *    If no job is set, means I should quit
 */
      if (!pp->jobf)
         break;
/*
 *    If there's work and I've got work to do
 */
      if (pp->chkin[rank] == 1)
      {
         pp->chkin[rank] = 2;
         ATL_tpool_dojob(pp, rank, 0);
      }
/*
 *    If no job, scope flag to see if I should sleep, or quit
 */
      else if (ATL_TPF_DIE(pp))
         break;
   }
   pp->chkin[rank] = 4;  /* tell master I'm completely done */
   #ifdef DEBUG
      fprintf(stderr, "%d: WORKER DONE\n", rank);
   #endif
   return(NULL);
}
#elif defined(ATL_POLLTPOOL) && !defined(ATL_OMP_THREADS)
/*
 * In this version of threadpool, worker threads restart master process once
 * job is done, but threads may still be awake and polling for work.  Since
 * the master process is not bound to a thread, this may cause contention
 * with the serial portion of the algorithm (if any), and so cause a slowdown.
 * However, it might be the right thing to do on systems where condition
 * variables are very expensive, or if the number of cores is very large,
 * or if polling does not interfere with the serial process for some reason.
 */
void *ATL_threadpool(void *vp)
{
   ATL_thread_t *tp=vp;
   ATL_tpool_t *pp = ATL_TP_PTR;
   const unsigned int rank=tp->rank;
   #ifdef ATL_PHI_SORTED
      const unsigned int nthr=pp->nthr, p4=nthr>>2, p4_2=p4+p4, p4_3=p4_2+p4;
   #endif
   void *tpmut = pp->tpmut;
   unsigned int jobID=0;

   #ifdef DEBUG
      fprintf(stderr, "%d: WORKER START P=%d\n", rank, pp->nthr);
   #endif
/*
 * Infinite loop over tasks.  Loop has 3 phases:
 * (1) Do work if commanded, including optional combine
 * (2) Poll for a period of time bounded by ATL_POLLTIME.  During this
 *     period the worker actively pools for work.  If a new job is requested
 *     for time expires, worker gets started on next job w/o th everhead
 *     of going to sleep & waking up.
 *     If ATL_POLLTIME <= 0, this phase is skipped.
 * (3) Sleep until master wakes up for new work
 * Have lock at top of loop.
 */
   ATL_mutex_lock(tpmut);  /* force wait until pp set up */
   while(1)
   {
      unsigned int vrank;
/*
 *    If I've got work to do
 */
      if (pp->jobf)
      {
         vrank = ATL_tpool_dojob(pp, rank, 0);
         ATL_mutex_lock(tpmut);
/*
 *       If I'm the last worker, signal master process we're done
 */
         if (vrank < pp->nworkers)
         {
            if (++pp->nwdone == pp->nworkers)
            {
               #ifdef DEBUG
                  fprintf(stderr, "%d,%d: SIGNAL MASTER\n", rank, vrank);
               #endif
               pp->NOWORK = pp->WORKDONE = 1;
               ATL_cond_signal(pp->mcond);
            }
         }
      }
/*
 *    If no job, scope flag to see if I should sleep, or quit
 */
      else
      {
         if (ATL_TPF_DIE(pp))
         {
            ATL_mutex_unlock(tpmut);
            return(NULL);
         }
      }
/*
 *    If we've finished so quickly that the job is still waiting on a
 *    mandatory worker, go back to job and pretend to be a new worker,
 *    as long as we are allowed to work on problem
 */
      #ifndef ATL_PHI_SORTED2
      #ifdef ATL_PHI_SORTED
         if (pp->wcnt < pp->nworkers && pp->nworkers != nthr &&
             pp->nworkers != p4_2 && pp->nworkers != p4_3)
      #else
         if (pp->wcnt < pp->nworkers)
      #endif
            continue;
      #endif
/*
 *    Sleep for next job; nsleep only examined during startup
 */
      pp->nsleep++;  /* don't care if this rolls over */
/*
 *    Poll until time elapses or a new job is available.
 */
      if (ATL_POLLTIME > 0.0 &&
         (jobID == pp->jobID || pp->wcnt >= pp->nworkers))
      {
         double t0, t1=0.0;
         ATL_mutex_unlock(tpmut);
         t0 = ATL_walltime();
         while (t1 < ATL_POLLTIME &&
                (jobID == pp->jobID || pp->wcnt >= pp->nworkers))
         {
            ATL_thread_yield();
            t1 = ATL_walltime() - t0;
         }
         ATL_mutex_lock(tpmut);
      }
/*
 *    If no new job available, go to sleep
 */
      while (pp->jobID == jobID || pp->wcnt >= pp->nworkers)
      {
         #ifdef DEBUG
            fprintf(stderr, "%d,%d: go to sleep\n", rank, vrank);
         #endif
         #ifdef ATL_PHI_SORTED
            ATL_cond_wait(pp->wcond, tpmut);
            if (rank < p4)
               ATL_cond_wait(pp->wcond, tpmut);
            else if (rank < p4_2)
               ATL_cond_wait(pp->wcond2, tpmut);
            else if (rank < p4_3)
               ATL_cond_wait(pp->wcond3, tpmut);
            else
               ATL_cond_wait(pp->wcond4, tpmut);
         #elif defined(ATL_PHI_SORTED2)
            if (rank >= p4_3)
               ATL_cond_wait(pp->wcond, tpmut)
            else if (rank < p4)
               ATL_cond_wait(pp->wconds[rank, tpmut]);
            else if (rank < p4_2)
               ATL_cond_wait(pp->wconds[rank-p4, tpmut]);
            else
               ATL_cond_wait(pp->wconds[rank-p4_2, tpmut]);
         #else
            ATL_cond_wait(pp->wcond, tpmut);
         #endif
      }  /* end sleep loop */
      jobID = pp->jobID;
   }     /* end loop over tasks */
/*
 * When a thread exits, decrement thread count
 */

   ATL_mutex_lock(tpmut);
   pp->nthr--;
   ATL_mutex_unlock(tpmut);
   #ifdef DEBUG
      fprintf(stderr, "%d: WORKER DONE\n", rank);
   #endif
   return(NULL);
}
#elif !defined(ATL_OMP_THREADS)
/*
 * In this version of threadpool, all worker threads go to sleep before
 * control is returned to the master process
 */
void *ATL_threadpool(void *vp)
{
   ATL_thread_t *tp=vp;
   ATL_tpool_t *pp = ATL_TP_PTR;
   const unsigned int rank=tp->rank;
   #if defined(ATL_PHI_SORTED) || defined(ATL_PHI_SORTED2)
      const unsigned int nthr=pp->nthr, p4=nthr>>2, p4_2=p4+p4, p4_3=p4_2+p4;
      const int icore = rank % (ATL_NTHREADS>>2);
   #endif
   void *tpmut = pp->tpmut;

   #ifdef DEBUG
      fprintf(stderr, "%d: WORKER START POLL=%d, P=%d\n", rank, DOPOLL, P);
   #endif
/*
 * Infinite loop over tasks.  Loop has 2 phases:
 * (1) Do work if commanded, including optional combine
 * (2) Wait for orders - can wait by polling or condition variable
 * Have lock at top of loop.
 */
   ATL_mutex_lock(tpmut);  /* force wait until pp set up */
   while(1)
   {
      unsigned int vrank;
/*
 *    If I've got work to do
 */
      if (pp->jobf)
      {
         vrank = ATL_tpool_dojob(pp, rank, 0);
         ATL_mutex_lock(tpmut);
/*
 *       Threads go to sleep or poll for work after each job
 *       If I'm the last guy, signal master process we're done
 */
         if (vrank < pp->nworkers)
         {
            if (++pp->nwdone == pp->nworkers)
            {
               pp->WORKDONE = 1;
               ATL_cond_signal(pp->mcond);
            }
         }
      }
/*
 *    If no work, scope flag to see if I should sleep, or quit
 */
      else
      {
         if (ATL_TPF_DIE(pp))
         {
            ATL_mutex_unlock(tpmut);
            return(NULL);
         }
      }
/*
 *    If the job is still waiting on workers, go back up and pretend to be
 *    a new worker as long as we're allowed to work on problem
 */
      #ifndef ATL_PHI_SORTED2
      #ifdef ATL_PHI_SORTED
         if (pp->wcnt < pp->nworkers && pp->nworkers != nthr &&
             pp->nworkers != p4_2 && pp->nworkers != p4_3)
      #else
         if (pp->wcnt < pp->nworkers)
      #endif
            continue;
      #endif
/*
 *    Sleep for next job; nsleep only examined during startup
 */
      pp->nsleep++;  /* don't care if this rolls over */
      do
      {
         #ifdef DEBUG
            fprintf(stderr, "%d: go to sleep\n", rank);
         #endif
         #ifdef ATL_PHI_SORTED
            if (rank < p4)
               ATL_cond_wait(pp->wcond, tpmut);
            else if (rank < p4+p4)
               ATL_cond_wait(pp->wcond2, tpmut);
            else if (rank < p4+p4+p4)
               ATL_cond_wait(pp->wcond3, tpmut);
            else
               ATL_cond_wait(pp->wcond4, tpmut);
         #elif defined(ATL_PHI_SORTED2)
            if (rank >= p4_3)
               ATL_cond_wait(pp->wcond, tpmut);
            else
               ATL_cond_wait(pp->wconds[icore], tpmut);
         #else
            ATL_cond_wait(pp->wcond, tpmut);
         #endif
      }
      while (pp->WORKDONE || pp->wcnt >= pp->nworkers);
   }  /* end loop over tasks */
/*
 * When a thread exits, decrement thread count
 */

   ATL_mutex_lock(tpmut);
   pp->nthr--;
   ATL_mutex_unlock(tpmut);
   #ifdef DEBUG
      fprintf(stderr, "%d: WORKER DONE\n", rank);
   #endif
   return(NULL);
}
#endif
