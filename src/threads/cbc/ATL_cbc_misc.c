#include <atlas_cbc.h>

void ATL_CBC_init(ATL_CBC_t *cbc, ATL_INT P, ATL_INT Q, ATL_CBC_TDATA *tdata)
{
   int i, k;
   assert(cbc);   /* check whether cbc is valid */
   assert(tdata); /* check whether tdata is valid */
   cbc->P = P;    /* number of rows in the thread grid */
   cbc->Q = Q;    /* number of columns in the thread grid */
   cbc->Nt = P*Q; /* total number of threads */
   cbc->tdata = tdata;  /* thread data that stores sync variables and ranks */

   cbc->Master = 0;     /* indicates who is the master for the grid */
   cbc->MallocOwner = 0;
   cbc->vp = malloc(sizeof(ATL_INT)*(P+Q)); /* memory for scoped masters */
   assert(cbc->vp);

   cbc->RowMasters = cbc->vp;       /* masters for each row */
   cbc->ColMasters = cbc->RowMasters + P; /* masters for each column */

   /* initialize the row masters to default */
   for (i=0, k=0; i<P; i++, k+=Q)
   {
      cbc->RowMasters[i] = k; /* row masters are threads of first column */
   }
   /* initialize the column masters to default */
   for (i=0; i<Q; i++)
   {
      cbc->ColMasters[i] = i; /* column masters are threads of first row */
   }
}

void ATL_CBC_destroy(ATL_CBC_t *cbc)
{
   if (cbc->vp)   /* if vp has not been freed yet, free the memory */
   {
      free((void *)cbc->vp);
      cbc->vp = NULL;   /* set vp to NULL to indicate that it's been freed. */
   }
}

void ATL_WaitForPartner(enum ATL_SYNC_SCOPE scope, ATL_CBC_t *cbc,
                        ATL_INT rankG, ATL_INT partner)
{
   while (cbc->tdata[rankG].NextMsg[scope]
            <= cbc->tdata[partner].NextMsg[scope]);
}

void ATL_PostNewMsg(enum ATL_SYNC_SCOPE scope, ATL_CBC_t *cbc,
                    ATL_INT rankG, ATL_INT msg)
{
   cbc->tdata[rankG].NextMsg[scope] = msg;
}

void ATL_WAIT(enum ATL_SYNC_SCOPE scope, ATL_CBC_t *cbc, ATL_INT rankG,
              ATL_INT waitOn, ATL_INT newMaster)
{
   if (waitOn == newMaster) /* No change in master. */
   {
      /* wait for master to finish. */
      while ((cbc->tdata[rankG]).Next[scope] !=
            (cbc->tdata[waitOn].Next[scope]))
            {
               ATL_tyield;
            }
   }
   else     /* Master will change. */
   {
      if (rankG == newMaster)  /* if I am the new master */
      {

         /* wait for old master to finish. */
         while ((cbc->tdata[rankG]).Next[scope] !=
            (cbc->tdata[waitOn]).Next[scope])
            {
               ATL_tyield;
            }

         /* Now update the master */
         if (scope == ATL_SYNC_GRID)   /* if whole grid */
         {
            cbc->Master = rankG;       /* update grid master */
         }
         else if (scope == ATL_SYNC_ROW)  /* if only a row */
         {
            int r = cbc->tdata[rankG].rankR; /* get my row rank */
            /* update the row master of my row to me */
            cbc->RowMasters[r] = rankG;
         }
         else if (scope == ATL_SYNC_COL)  /* if only a column */
         {
            int c = cbc->tdata[rankG].rankC; /* get my col rank */
            /* update the row master of my row to me */
            cbc->ColMasters[c] = rankG;
         }
      }
      else /* I'm not the new master,wait for the new master to take control. */
      {
         if (scope == ATL_SYNC_GRID)   /* if whole grid */
         {
            /* wait for grid master to change */
            while (cbc->Master != newMaster)
            {
               ATL_tyield;
            }
         }
         else if (scope == ATL_SYNC_ROW)  /* if only a row */
         {
            int r = cbc->tdata[rankG].rankR; /* get my row rank */
            /* wait for row master to change */
            while (cbc->RowMasters[r] != newMaster)
            {
               ATL_tyield;
            }
         }
         else if (scope == ATL_SYNC_COL)  /* if only a column */
         {
            int c = cbc->tdata[rankG].rankC; /* get my col rank */
            /* wait for colum master to change */
            while (cbc->ColMasters[c] != newMaster)
            {
               ATL_tyield;
            }
         }
      }
   }
}


