#include "atlas_misc.h"
#include "atlas_threads.h"
/*
 * GlobalAtomicCountDowns are same as GlobalAtomicCount, except that they
 * guarantee that their last non-zero value is 1 (they don't guarantee
 * a sequential count down, despite the name).  They are used over Counts
 * when it is vital that the someone atomically knows they are the last
 * non-zero return (i.e. the last worker frees resources, etc).
 */
void *ATL_SetGlobalAtomicCountDown
(
   int P,               /* # of Local counters to use to make global ctr */
   int cnt              /* total count to start global count at */
)
{
   void **va;
   va = malloc(2*sizeof(void*));
   ATL_assert(va);
   if (cnt > 2 && P > 1)
   {
      va[0] = ATL_SetGlobalAtomicCount(P, cnt-1, 0);
      va[1] = ATL_SetAtomicCount(1);
   }
   else
   {
      va[0] = NULL;
      va[1] = ATL_SetAtomicCount(cnt);
   }
   return(va);
}

