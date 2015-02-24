#include "atlas_misc.h"
#include "atlas_threads.h"
#ifndef ATL_WINTHREADS
   #include <sched.h>
#endif
void ATL_thread_yield(void)
{
   #ifdef ATL_WINTHREADS
      Sleep(0);
   #elif defined(ATL_ARCH_XeonPHI)
       __asm__ __volatile__ ("delay %0" :: "r"(64) : "0" );
   #else
      sched_yield();
   #endif
}
