#include "atlas_misc.h"
#include "atlas_threads.h"
/*
 * condition variable funcs are used in some threadpool implementations, but
 * parallel section sync handles this for OpenMP,
 * On Windows we do not presently support (must use locks or FULLPOLL),
 * so simply issue runtime assert(0) for these implementations so we can
 * find bad calls on Windows or OpenMP.
 */
#if defined(ATL_OMP_THREADS) || defined(ATL_WINTHREADS)
   #define ATL_DIE 1
#endif
void *ATL_cond_init(void)
{
   #ifdef ATL_DIE
      ATL_assert(0);
      return(NULL);
   #else
      void *vp;
      vp = malloc(sizeof(pthread_cond_t));
      ATL_assert(vp);
      ATL_assert(!pthread_cond_init(vp, NULL));
      return(vp);
   #endif
}
