#include "atlas_misc.h"
#include "atlas_threads.h"

void ATL_goparallel
/*
 * This function is used when you pass a single opstruct to all threads;
 * In this case, we stash opstruct in launchstruct's vp, and then use the
 * rank array as opstruct during the spawn.  Therefore, these routines
 * should expect to get their problem def from ls.vp, and their rank from
 * the second argument.  The DoWork function is the function that should
 * be called from each thread to do the parallel work.  This function should
 * look like:
 * void DoWork_example(ATL_LAUNCHSTRUCT_t *lp, void *vp)
 * {
 *    ATL_thread_t *tp = vp;
 *    const int myrank = tp->rank;
 *    my_prob_def_t *pd = lp->opstruct;
 *    ... do work based on info in struct pointed to by lp->opstruct ...
 * }
 * Your DoWork should perform any needed combine before finishing execution,
 * and any return values can be passed in the problem definition structure
 * that you define.
 */
(
   const unsigned int P, /* # of cores to use */
   void *DoWork,         /* func ptr to work function */
   void *opstruct,       /* structure giving tasks to threads */
   void *DoComb          /* function to combine two opstructs */
)
{
   ATL_thread_t *tp;
   int *chkin;
   void *vp, *lc;
   int i;
   ATL_LAUNCHSTRUCT_t ls;

   ls.OpStructIsInit = NULL;
   ls.DoWork = DoWork;
   ls.DoComb = DoComb;
   ls.opstruct = opstruct;
   if (DoComb)
      ATL_goParallel(P, ATL_oldjobwrap, ATL_oldcombwrap, &ls, NULL);
   else
      ATL_goParallel(P, ATL_oldjobwrap, NULL, &ls, NULL);
}
