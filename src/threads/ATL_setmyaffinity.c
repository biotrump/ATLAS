#include "atlas_taffinity.h"  /* include this file first! */
/*
 * Cases handled before thread creation / start running:
 *    ATL_WINTHREADS, ATL_PAFF_SETAFFNP, ATL_PAFF_SETPROCNP
 * Cases for post-creation binding done by spawner:
 *    ATL_PAFF_PBIND, ATL_PAFF_BINDP, ATL_PAFF_CPUSET
 * all of these cases provide self-setting as well.
 * Cases that work only for self-setting:
 *   ATL_PAFF_SCHED, ATL_PAFF_PLPA
 */

#ifdef ATL_WINTHREADS
   #include <windows.h>
#elif defined(ATL_NOAFFINITY)
   #ifndef ATL_OMP_THREADS
      #include <pthread.h>
   #endif
#elif defined(ATL_SPAFF_PLPA)
   #include <pthread.h>
   #include <plpa.h>
#elif defined(ATL_SPAFF_PBIND)
   #include <pthread.h>
   #include <sys/types.h>
   #include <sys/processor.h>
   #include <sys/procset.h>
#elif defined(ATL_SPAFF_SCHED)
   #define _GNU_SOURCE 1 /* what manpage says you need to get CPU_SET */
   #define __USE_GNU   1 /* what works on linuxes that I've seen */
   #include <sched.h>    /* must include this before pthreads */
   #include <pthread.h>
#elif defined(ATL_SPAFF_RUNON)
   #include <pthread.h>
#elif defined(ATL_SPAFF_BINDP)
   #include <pthread.h>
   #include <sys/thread.h>      /* thread_self header */
   #include <sys/processor.h>   /* bindprocessor header */
#elif defined(ATL_SPAFF_CPUSET)
   #include <pthread.h>
   #include <sys/param.h>
   #include <sys/cpuset.h>
#endif
#include "atlas_misc.h"
#include "atlas_threads.h"
int ATL_setmyaffinity(const int rank)
/*
 * Attempts to sets the affinity of an already-running thread.  The
 * aff_set flag is set to true whether we succeed or not (no point in
 * trying multiple times).
 * RETURNS: 0 on success, non-zero error code on error
 */
{
   #if !defined(ATL_AFF_NUMID)
      const int bindID=rank;
   #elif defined(ATL_RANK_IS_PROCESSORID)
      const int bindID = rank % ATL_AFF_NUMID;
   #else
      const int bindID = ATL_affinityIDs[rank%ATL_AFF_NUMID];
   #endif
#ifdef ATL_WINTHREADS  /* Windows */
      ATL_assert(SetThreadAffinityMask(GetCurrentThreadId(),
                                       (((long long)1)<<bindID)));
#elif defined(ATL_SPAFF_PLPA)  /* affinity wrapper package */
   plpa_cpu_set_t cpuset;
   PLPA_CPU_ZERO(&cpuset);
   PLPA_CPU_SET(bindID, &cpuset);
   return(plpa_sched_setaffinity((pid_t)0, sizeof(cpuset), &cpuset));
#elif defined(ATL_SPAFF_PBIND)
   return(processor_bind(P_LWPID, P_MYID, bindID, NULL));
#elif defined(ATL_SPAFF_SCHED)  /* linux */
   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(bindID, &cpuset);
   return(sched_setaffinity(0, sizeof(cpuset), &cpuset));
#elif defined (ATL_SPAFF_RUNON)  /* IRIX: this will not work for process! */
   return(pthread_setrunon_np(bindID));
#elif defined(ATL_SPAFF_BINDP)  /* AIX */
   return(bindprocessor(BINDTHREAD, thread_self(), bindID));
#elif defined(ATL_SPAFF_CPUSET)  /* untried FreeBSD code */
   cpuset_t mycpuset;
   CPU_ZERO(&mycpuset);         /* no manpage, so guess works like linux */
   CPU_SET(bindID, &mycpuset);
   if (me->affID >= 0)
      return(0);
   me->affID = bindID;
   return(cpuset_setaffinity(CPU_LEVEL_WHICH, CPU_WHICH_TID, -1,
                             sizeof(mycpuset), &mycpuset));
#elif defined(ATL_NOAFFINITY)
   return(0);   /* no-affinity override says not setting is not error */
#else
   return(1);  /* Don't know how to set affinity, return error */
#endif
   return(0);
}
