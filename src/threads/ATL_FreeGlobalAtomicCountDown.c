#include "atlas_misc.h"
#include "atlas_threads.h"
/*
 * GlobalAtomicCountDowns are same as GlobalAtomicCount, except that they
 * guarantee that their last non-zero value is 1 (they don't guarantee
 * a sequential count down, despite the name).  They are used over Counts
 * when it is vital that the someone atomically knows they are the last
 * non-zero return (i.e. the last worker frees resources, etc).
 */
void ATL_FreeGlobalAtomicCountDown(void *vp)
{
   void **va=vp;
   if (va[0])
      ATL_FreeGlobalAtomicCount(va[0]);
   ATL_FreeAtomicCount(va[1]);
   free(vp);
}
