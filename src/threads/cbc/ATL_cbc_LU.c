#include "atlas_cbc.h"
void ATL_CBC_LUschedWork(enum ATL_SYNC_SCOPE scope, ATL_CBC_t *cbc,
                         ATL_INT rankG, ATL_INT *next_val, ATL_INT atSwitch)
{
#if 1
   int m_val = *next_val + 1;
   int i, m_val2, _needToWait;
   int rid = cbc->tdata[rankG].rankR;
   int cid = cbc->tdata[rankG].rankC;
   //_T("T%d: Trying barrier on %d\n", rankG, m_val);
   cbc->tdata[rankG].NextMsg[scope] = m_val;
   if (scope == ATL_SYNC_COL)
   {
      for (i=cid; i<cbc->Nt; i+=cbc->Q)
      {
         if (i != rankG)
         {
            while(cbc->tdata[i].NextMsg[scope] <= 0)
            {
               ATL_tyield;
            }
            cbc->tdata[rankG].NextMsg[scope] =
               Mmin(cbc->tdata[rankG].NextMsg[scope],
                  cbc->tdata[i].NextMsg[scope]);
         }
      }
      if (m_val > cbc->tdata[rankG].NextMsg[scope]) /* i was in wrong column */
      {
         *next_val = cbc->tdata[rankG].NextMsg[scope] - 1;
         //_T("T%d: Exiting wrong barrier on %d\n", rankG, m_val);
         return; /* return and come back with correct column */
      }
      else
      {
         ATL_CBC_barrier(scope, cbc, rankG);
      }
   }
   //_T("T%d: Exiting barrier on %d\n", rankG, m_val);
   cbc->tdata[rankG].NextMsg[scope] = -1;
   ATL_CBC_barrier(scope, cbc, rankG);
#else
   int m_val = *next_val + m_val + 1;
   int i,  m_val2, _needToWait;
   int rid = cbc->tdata[rankG].rankR;
   int cid = cbc->tdata[rankG].rankC;
   while (1)
   {
      cbc->tdata[rankG].NextMsg[scope] = m_val;
      _needToWait = 0;
      if (scope == ATL_SYNC_COL)
      {
         for (i=cid; i<cbc->Nt; i+=cbc->Q)
         {
            if (i != rankG)
            {
               while(cbc->tdata[i].NextMsg[scope] <= 0);
               if (cbc->tdata[rankG].NextMsg[scope] >
                     cbc->tdata[i].NextMsg[scope])
               {
                  /* i am in wrong column */
                  cbc->tdata[rankG].NextMsg[scope] =
                     cbc->tdata[i].NextMsg[scope];
               }
               else if (cbc->tdata[rankG].NextMsg[scope] <
                     cbc->tdata[i].NextMsg[scope])
               {
                  _needToWait = 1; /* someone else was in wrong column */
               }
            }
         }
         m_val2 = cbc->tdata[rankG].NextMsg[scope];
         ATL_CBC_barrier(scope, cbc, rankG);
         cbc->tdata[rankG].NextMsg[scope] = -1;
         ATL_CBC_barrier(scope, cbc, rankG);
         if (m_val > m_val2) /* i was in wrong column */
         {
            *next_val = m_val2 - 1;
            //_T("T%d: Exiting wrong barrier on %d\n", rankG, m_val);
            return; /* return and come back with correct column */
         }
         else if (_needToWait)
            /* someone else was in wrong panel, need to wait for them */
         {
            continue;
         }
         break;
      }
   }
#endif
}



void ATL_iTMaxForLU(enum ATL_SYNC_SCOPE scope, ATL_CBC_t *cbc, ATL_INT rankG,
                    ATL_INT master, ATL_INT newMaster, const TYPE in_data,
                    const int in_indx, TYPE *out_data, int *out_indx, int *ownr)
{
   ATL_TDATA(rankG).mdata = in_data;
   ATL_TDATA(rankG).mindx = in_indx;
   ATL_TDATA(rankG).mOwner = rankG;

   if (scope == ATL_SYNC_GRID)   /* Combine is for whole grid */
   {
      int d = 1, partner, vrankG;        /* d is the partner distance */
      #ifdef USE_ASSERT
         ATL_assert(master == cbc->Master);
      #endif
      vrankG = rankG - master + cbc->Nt;   /* find the virtual rank : shift */
      vrankG = vrankG - (vrankG / cbc->Nt) * cbc->Nt; /* vrankG %= Nt */
      while ((vrankG & d) == 0 &&   /* if I need to wait on someone and that */
            (vrankG + d) < cbc->Nt) /* someone exists, enter this loop */
      {
         partner = vrankG + d + master;      /* get partner's real rank */
         partner = partner - (partner / cbc->Nt) * cbc->Nt; /* partner %= Nt */
         /* wait for partner to finish */
         while (ATL_NEXT(cbc, rankG, scope) == ATL_NEXT(cbc, partner, scope));

         /* do the combine (max loc for this case) */
         if (Mabs(ATL_TDATA(rankG).mdata) < Mabs(ATL_TDATA(partner).mdata))
         {
            ATL_TDATA(rankG).mdata = ATL_TDATA(partner).mdata;
            ATL_TDATA(rankG).mindx = ATL_TDATA(partner).mindx;
            ATL_TDATA(rankG).mOwner = ATL_TDATA(partner).mOwner;
         }

         d <<= 1;                   /* increase partner distance */
      } /* end of while loop of binary combine */
   } /* if statement for grid scope is done */

   else if (scope == ATL_SYNC_ROW) /* for row scope */
   {
      int d = 1, partner, vrankC;
      #ifdef USE_ASSERT
      ATL_assert(master == cbc->RowMasters[ATL_TDATA(rankG).rankR]);
      ATL_assert(ATL_ROW_RANK(cbc, rankG) == ATL_ROW_RANK(cbc, master));
      ATL_assert(ATL_ROW_RANK(cbc, newMaster) == ATL_ROW_RANK(cbc, master));
      #endif
      vrankC = ATL_COL_RANK(cbc, rankG) - ATL_COL_RANK(cbc, master)
               + cbc->Q;
      vrankC = vrankC - (vrankC / cbc->Q) * cbc->Q;
      while ((vrankC &d) == 0 &&
            (vrankC + d) < cbc->Q)
      {
         partner = vrankC + d + ATL_COL_RANK(cbc, master); /* real column */
         partner = partner - (partner / cbc->Q) * cbc->Q;
         partner = ATL_GRID_RANK(cbc,
               ATL_ROW_RANK(cbc, master), partner); /*grid rank*/

         /* wait for partner to finish */
         while (ATL_NEXT(cbc, rankG, scope) == ATL_NEXT(cbc, partner, scope));

         /* do the combine (max loc for this case) */
         if (Mabs(ATL_TDATA(rankG).mdata) < Mabs(ATL_TDATA(partner).mdata))
         {
            ATL_TDATA(rankG).mdata = ATL_TDATA(partner).mdata;
            ATL_TDATA(rankG).mindx = ATL_TDATA(partner).mindx;
            ATL_TDATA(rankG).mOwner = ATL_TDATA(partner).mOwner;
         }

         d <<= 1;                   /* increase partner distance */
      } /* end of while loop of binary combine */
   } /* if statement for row scope is done */

   else  /* for column scope */
   {
      int d = 1, partner, vrankR;
      #ifdef USE_ASSERT
         ATL_assert(master == cbc->ColMasters[ATL_TDATA(rankG).rankC]);
         ATL_assert(ATL_COL_RANK(cbc, rankG) == ATL_COL_RANK(cbc, master));
         ATL_assert(ATL_COL_RANK(cbc, newMaster) == ATL_COL_RANK(cbc, master));
      #endif
      vrankR = ATL_ROW_RANK(cbc, rankG) - ATL_ROW_RANK(cbc, master)
               + cbc->P;
      vrankR = vrankR - (vrankR / cbc->P) * cbc->P;
      while ((vrankR &d) == 0 &&
            (vrankR + d) < cbc->P)
      {
         partner = vrankR + d + ATL_ROW_RANK(cbc, master); /* real row */
         partner = partner - (partner / cbc->P) * cbc->P;
         partner = ATL_GRID_RANK(cbc, partner,
               ATL_COL_RANK(cbc, master)); /*grid rank*/

         /* wait for partner to finish */
         while (ATL_NEXT(cbc, rankG, scope) == ATL_NEXT(cbc, partner, scope));

         /* do the combine (max loc for this case) */
         if (Mabs(ATL_TDATA(rankG).mdata) < Mabs(ATL_TDATA(partner).mdata))
         {
            ATL_TDATA(rankG).mdata = ATL_TDATA(partner).mdata;
            ATL_TDATA(rankG).mindx = ATL_TDATA(partner).mindx;
            ATL_TDATA(rankG).mOwner = ATL_TDATA(partner).mOwner;
         }

         d <<= 1;                   /* increase partner distance */
      } /* end of while loop of binary combine */
   } /* if statement for colum scope is done */

   /* All done with binary combine, now update state */
   if (rankG == master) /* if I am master */
   {
#if 1    /* leave on all, set to 0 if leave on only the master */

      int i, j, k, first, limit;
      k = 1;                     /* for grid and column scope, increment 1 */
      first = 0;                 /* for grid, start from 0 */
      limit = cbc->Nt;           /* how many threads to check, for grid, all */

      if (scope == ATL_SYNC_COL) /* if column scope */
      {
         k = cbc->Q;             /* increment is number of columns */
         limit = cbc->P;         /* no. of threads to check is no. of rows */
         first = (cbc->tdata[rankG]).rankC;  /* first is just column rank */
      }
      else if (scope == ATL_SYNC_ROW) /* if row scope */
      {
         limit = cbc->Q;         /* no. of threads to check is no. of cols */

         /* first is the first of that row */
         first = (cbc->tdata[rankG]).rankR * cbc->Q;
      }

      /* now copy my result to everyone else */
      for (i=0, j=first; i<limit; i++, j+=k)
      {
         ATL_TDATA(j).mdata = ATL_TDATA(rankG).mdata;
         ATL_TDATA(j).mindx = ATL_TDATA(rankG).mindx;
         ATL_TDATA(j).mOwner = ATL_TDATA(rankG).mOwner;
      }
#endif

      /* POST so that the everyone knows that I am done. */
      ATL_POST(scope, cbc, rankG);
      if (master != newMaster) /* I am not the new master */
         /* wait for new master*/
         ATL_WAIT(scope, cbc, rankG, master, newMaster);
   }
   else /* I am not old master */
   {
      /* done my part of combine */
      ATL_POST(scope, cbc, rankG);  /* let others know that I am done */
      ATL_WAIT(scope, cbc, rankG, master, newMaster); /* wait for newMaster */
   }
   *out_data = ATL_TDATA(rankG).mdata;
   *out_indx = ATL_TDATA(rankG).mindx;
   *ownr = ATL_TDATA(rankG).mOwner;
}

void ATL_CBC_Malloc_Lock(ATL_CBC_t *cbc, ATL_INT rankG)
{
   while (cbc->MallocOwner != rankG);
}

void ATL_CBC_Malloc_UnLock(ATL_CBC_t *cbc, ATL_INT rankG)
{
   assert(cbc->MallocOwner == rankG);
   cbc->MallocOwner = (cbc->MallocOwner+1)%cbc->Nt;
}
