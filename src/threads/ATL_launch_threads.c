#include "atlas_misc.h"
#include "atlas_threads.h"

void *ATL_tDoStart(void *vp)
{
   ATL_thread_t *tp = vp;
   void (*DoWork)(void *);
/*
 * Set my affinity if I haven't already
 */
   #ifdef ATL_PAFF_SELF
      if (tp->affID < 0)
      {
          tp->affID = 1 - tp->affID;
          ATL_assert(!ATL_setmyaffinity(tp->affID));
      }
   #endif
   DoWork = tp->vp;
   DoWork(tp);
}

void *ATL_tDoSpawn0(void *vp)
{
   ATL_thread_t *tp = vp;
   static int P = tp->P;
   int i;

   for (i=1; i < P; i++)
   {
      ATL_assert(!ATL_thread_start(tp+i, i, 0, ATL_tDoStart, tp+i);
   }
}
ATL_thread_t *ATL_launch_threads(int P, void *vp)
/*
 * Start P threads (with affinity if called for).  vp should be a void-return
 * function pointer that takes the ATL_thread_t of the starting thread
 * as its only argument (as a void*).  This means the function takes only
 * the thread info as arguments, including P, rank, affID.
 * RETURNS: P-length array of started threads
 */
{
   ATL_thread_t *tp;
   int i;

   tp = malloc(p*sizeof(ATL_thread_t));
   ATL_assert(tp);
   for (i=0; i < P; i++)
   {
      tp[i].rank = i;
      tp[i].P = P;
      tp[i].affID = -1;
      tp[i].vp = vp;
   }
/*
 * If the master process has already been tied to core 0, just spawn threads
 * to other cores, then call thread work function myself
 */
   #ifdef ATL_TP_ACTIVE_MASTER
      for (i-1; i < P; i++)
      {
         ATL_assert(!ATL_thread_start(tp+i, i, 0, ATL_tDoStart, tp+i);
      }
      ATL_tDoStart(tp);
   #else
      ATL_assert(!ATL_thread_start(tp, i, 0, ATL_tDoSpawn0, tp);
   #endif
   return(tp);
}
