/*
 * This file implements barrier primitives using CBC
 * (see atlas_cbc.h for details).
 *
 * Originally it was done using a counter but to avoid overflow we used
 * boolean state here. There was a race condition when using == and != rather
 * than using >= and > . To avoid this, Anthony Castaldo suggested the idea
 * of changing the masters.
 *
 * NOTE: We also tried to use the xor of rank aka virtual rank to implement
 * the barrier. But this has the similar issue as before but this time the
 * problem occurs when the master is not changing.
 */

#include "atlas_cbc.h"
#include <string.h>


#ifndef Mabs
   #define Mabs(a) ((a) < 0 ? -(a) : (a))
#endif

#if 0
   #define USE_ASSERT
#endif

void ATL_CBC_barrier_init(ATL_CBC_t *cbc, int pt, int qt)
{
   const int nt = pt * qt;
   int i, j, k;

   ATL_CBC_TDATA* tdata = malloc(sizeof(ATL_CBC_TDATA)*nt);
   for (i=0, k=0; i<pt; i++)
   {
      for (j=0; j<qt; j++)
      {
         tdata[k].Next[0] = 0;
         tdata[k].Next[1] = 0;
         tdata[k].Next[2] = 0;
         tdata[k].NextMsg[0] = 0;
         tdata[k].NextMsg[1] = 0;
         tdata[k].NextMsg[2] = 0;
         tdata[k].rankR = i;
         tdata[k].rankC = j;
         tdata[k].rankG = k++;
      }
   }
   ATL_CBC_init(cbc, pt, qt, tdata);
}

void ATL_CBC_barrier_destroy(ATL_CBC_t *cbc)
{
   if (cbc && cbc->tdata)
   {
      free(cbc->tdata);
      cbc->tdata = NULL;
   }
   ATL_CBC_destroy(cbc);
}

void ATL_CBC_barrier_internal(enum ATL_SYNC_SCOPE scope, ATL_CBC_t *cbc,
                              ATL_INT rankG, ATL_INT waitOn, ATL_INT newMaster)
{
   if (rankG == waitOn)
   {
      int i, j, k, first, limit;
      k = 1;                     /* for grid and column scope, increment 1 */
      first = 0;                 /* for grid, start from 0 */
      limit = cbc->Nt;           /* how many threads to check, for grid, all */

      if (scope == ATL_SYNC_COL) /* if column scope */
      {
         #ifdef USE_ASSERT
            ATL_assert(waitOn == cbc->ColMasters[ATL_TDATA(rankG).rankC]);
         #endif
         k = cbc->Q;             /* increment is number of columns */
         limit = cbc->P;         /* no. of threads to check is no. of rows */
         first = (cbc->tdata[rankG]).rankC;  /* first is just column rank */
      }
      else if (scope == ATL_SYNC_ROW) /* if row scope */
      {
         #ifdef USE_ASSERT
            ATL_assert(waitOn == cbc->RowMasters[ATL_TDATA(rankG).rankR]);
         #endif
         limit = cbc->Q;         /* no. of threads to check is no. of cols */

         /* first is the first of that row */
         first = (cbc->tdata[rankG]).rankR * cbc->Q;
      }
      else
      {
         #ifdef USE_ASSERT
            ATL_assert(waitOn == cbc->Master);
         #endif
      }

      /* now wait for all required threads to finish */
      for (i=0, j=first; i<limit; i++, j+=k)
      {
         if (j != rankG)      /* if it's not me, wait for it */
         {
            while ((cbc->tdata[rankG]).Next[scope] ==
                     (cbc->tdata[j]).Next[scope])
            {
               ATL_tyield;
            }
         }
      }
      /* All threads done, now I can post. */
      ATL_POST(scope, cbc, rankG);

      /* if I am not the new Master, wait for new master to take over */
      if (waitOn != newMaster)
         ATL_WAIT(scope, cbc, rankG, waitOn, newMaster);
   }
   else  /* I am not the master */
   {
      ATL_POST(scope, cbc, rankG);  /* post my update */

      /* wait for old to post update or new master to take over */
      ATL_WAIT(scope, cbc, rankG, waitOn, newMaster);
   }
}

void ATL_CBC_barrier(enum ATL_SYNC_SCOPE scope, ATL_CBC_t *cbc, ATL_INT rankG)
{
   int m;
   //_T("T%d: entering barrier ...\n", rankG);
   switch (scope)
   {
      case ATL_SYNC_GRID:
         ATL_CBC_barrier_internal(scope, cbc, rankG, cbc->Master, cbc->Master);
         break;
      case ATL_SYNC_ROW:
         m = cbc->tdata[rankG].rankR; /* get my col rank */
         m = cbc->RowMasters[m];
         ATL_CBC_barrier_internal(scope, cbc, rankG, m, m);
         break;
      case ATL_SYNC_COL:
         m = cbc->tdata[rankG].rankC; /* get my col rank */
         m = cbc->ColMasters[m];
         ATL_CBC_barrier_internal(scope, cbc, rankG, m, m);
         break;
   }
   //_T("T%d: exiting barrier ...\n", rankG);
}
void ATL_CBC_gbarrier(ATL_CBC_t *cbc, int rankG)
{
   ATL_CBC_barrier(ATL_SYNC_GRID, cbc, rankG);
}

