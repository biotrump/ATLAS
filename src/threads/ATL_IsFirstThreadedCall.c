#include "atlas_threads.h"
#ifdef ATL_OMP_THREADS
   #include "atlas_misc.h"
#endif
/*
 * A kludge for initializing thread pool related stuff on first call
 * in portable fashion.
 * OpenMP just uses the OpenMP thread pool, so should not call here.
 * pthreads uses standard initializer.
 * Otherwise, I assume Windows, where DecAtomicCtr should work.
 * Note that this assumption is wrong for Windows on ARM (god help us).
 */
#if ATL_USE_POSIXTHREADS
   static pthread_mutex_t initlock=PTHREAD_MUTEX_INITIALIZER;
   static int atlinit=1;
#else   /* non-pthreads is Windows, where AtomicCtrs work */
   volatile static int ilck[65] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1};
#endif
int ATL_IsFirstThreadedCall(void)
{
#ifdef ATL_OMP_THREADS
   ATL_assert(0);  /* OpenMP should never call this! */
#elif ATL_USE_POSIXTHREADS
   int iret;
   pthread_mutex_lock(&initlock);
   iret = atlinit;
   atlinit = 0;
   pthread_mutex_unlock(&initlock);
   return(iret);
#else   /* non-pthreads is Windows, where AtomicCtrs work */
   return(ATL_DecAtomicCtr((void*)ilck);
#endif
}
/*
 * This function not thread safe, so user must make sure only called once
 */
void ATL_ResetThreadPoolInit(void)
{
#if ATL_USE_POSIXTHREADS
   pthread_mutex_lock(&initlock);
   atlinit = 1;
   pthread_mutex_unlock(&initlock);
#else
   ATL_ResetAtomicCount((void*)ilck, 1);
#endif
}
