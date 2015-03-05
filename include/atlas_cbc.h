#ifndef ATLAS_CBC_H
#define ATLAS_CBC_H
/*
 * This file prototypes cache-based communication (CBC) routines.  CBC
 * exploits cache coherence protocols for synchronization and communication
 * that runs at the speed of hardware, rather than the speed of software
 * that you potentially get when calling OS-provided functions like mutex, etc.
 * CBC requires strongly-ordred cache choerence to work correctly.  In
 * a strongly ordered cache, if a given core writes memory two memory
 * locations A and B in that sequence, the coherence protocal guarantees
 * that other cores only "see" the change in B *after* they can "see"
 * the change in A.  In weakly ordered caches, the core may "see" the
 * changes in any order (so seeing B change does not guarantee you won't
 * get a stale value if you then read A).  In general, x86 architectures
 * have strongly ordered caches, while non-x86 (eg., ARM and PowerPC) have
 * weakly-ordered caches.  Therefore, most CBC routines should not be used
 * on non-x86 platforms unless they are combined with memory fences of some
 * sort.  On x86 at least, CBC is critical to driving down communication costs
 * on extreme parallel architectures like the Xeon PHI.
 */
#if 0
   #include "atlas_threads.h"
   #include "atlas_lapack.h"
#endif
#include <assert.h>
#include <malloc.h>
#include <stdio.h>
#include <unistd.h>
#include "atlas_misc.h"

#if 1
#ifndef TYPE
   #define ATL_INT int
   #define TYPE double
#endif
#endif

#if 1 && defined(ATL_ARCH_XeonPHI)
   #define ATL_tyield __asm__ __volatile__ ("delay %0" :: "r"(64) : "0" )
#elif 1 && defined(USE_YIELD)
   #define ATL_tyield pthread_yield()
#else
   #define ATL_tyield
#endif

enum ATL_SYNC_SCOPE { ATL_SYNC_GRID=0, ATL_SYNC_ROW=1, ATL_SYNC_COL=2 };
enum ATL_CBC_OP { ATL_CBC_SUM, ATL_CBC_MAX, ATL_CBC_MIN, ATL_CBC_MAXLOC };
enum ATL_CBC_TYPE { ATL_CBC_INT=0,
                    ATL_CBC_FLOAT=1,
                    ATL_CBC_DOUBLE=2,
                    ATL_CBC_FLOAT_INT=3,
                    ATL_CBC_DOUBLE_INT=4
                  };

static int ATL_CBC_TYPE_SIZE[5] = { sizeof(int), sizeof(float), sizeof(double),
                              sizeof(int)+sizeof(float),
                              sizeof(int)+sizeof(double)
                           };

typedef struct
{
   unsigned int rankG;      /* global rank */
   unsigned int rankR;      /* row rank */
   unsigned int rankC;      /* column rank */
   //int ierr;
   void *workspace;  /* ptr to array of pointers */
   int *ldws;
   void *progs;      /* ptr to array of pointers for progress */

   /* Following data shouldn't be accessed by user code, only CBC */
   volatile int Next[3];    /* 3 scopes, grid, row, column */
   volatile int NextMsg[3];   /* 3 scopes, these will be incremented */

#if 0
   enum ATL_CBC_TYPE ctype;   /* data type for combine */
   void *cdata_in;            /* input data for combine */
   ATL_INT clen_in;           /* data length for combine */
   void *cdata_out;           /* output data for combine */
   ATL_INT clen_out;          /* output data length */
#endif

   volatile TYPE mdata;
   volatile int mindx;
   volatile int mOwner;

   void *cta_args;            /* context-aware args */
   volatile int cta_state[1]; /* context-aware state */
   void *cbc_bar;             /* local CBC variable */
   //void *pt_bar;            /* pthread_barrier for now */

   char space[128
      - sizeof(void*)            /* cta_args */
      - sizeof(volatile int)*1   /* cta_state */
      - sizeof(unsigned int) /* syncCount */
      - 3*sizeof(unsigned int)
      //- sizeof(int)
      /*
      - sizeof(ATL_INT)*2
      - sizeof(void *)*2
      - sizeof(enum ATL_CBC_TYPE)
      */
      - sizeof(void*)
      - sizeof(int*)
      - sizeof(void*)
      - 3*sizeof(volatile int)
      - 3*sizeof(volatile int)
      - sizeof(TYPE)
      - 2*sizeof(int)];/* Available space...to fill cache line */
} ATL_CBC_TDATA;

typedef struct
{
   float fdata;   /* data that needs to be compared */
   ATL_INT index; /* index for that data */
} ATL_CBC_FLOAT_INT_t;

typedef struct
{
   double fdata;   /* data that needs to be compared */
   ATL_INT index; /* index for that data */
} ATL_CBC_DOUBLE_INT_t;

typedef struct
{
   ATL_INT P;
   ATL_INT Q;
   ATL_INT Nt;
   ATL_CBC_TDATA *tdata;
   volatile ATL_INT Master;
   volatile ATL_INT *RowMasters;
   volatile ATL_INT *ColMasters;
   void *vp;
   volatile ATL_INT MallocOwner;
} ATL_CBC_t;

void ATL_CBC_barrier_init(ATL_CBC_t*, int, int);
void ATL_CBC_barrier_destroy(ATL_CBC_t*);
void ATL_CBC_gbarrier(ATL_CBC_t*, int);

void ATL_CBC_init(ATL_CBC_t*, ATL_INT, ATL_INT, ATL_CBC_TDATA*);
void ATL_CBC_destroy(ATL_CBC_t*);
void ATL_WAIT(enum ATL_SYNC_SCOPE, ATL_CBC_t*, ATL_INT, ATL_INT, ATL_INT);
void ATL_CBC_barrier_internal(enum ATL_SYNC_SCOPE, ATL_CBC_t*, ATL_INT,
                     ATL_INT, ATL_INT);
void ATL_CBC_barrier(enum ATL_SYNC_SCOPE, ATL_CBC_t*, ATL_INT);
void ATL_CBC_LUschedWork(enum ATL_SYNC_SCOPE, ATL_CBC_t*,
      ATL_INT,ATL_INT*,ATL_INT);

void ATL_Send(ATL_CBC_t*, ATL_INT, ATL_INT, void*, ATL_INT);
void ATL_Recv(ATL_CBC_t*, ATL_INT, ATL_INT, void*, ATL_INT);

void ATL_Combine(enum ATL_SYNC_SCOPE, ATL_CBC_t*, ATL_INT, ATL_INT, ATL_INT,
                 enum ATL_CBC_TYPE, enum ATL_CBC_OP, void*, ATL_INT, void*);

void ATL_iTMaxForLU(enum ATL_SYNC_SCOPE, ATL_CBC_t*, ATL_INT, ATL_INT, ATL_INT,
                    const TYPE, const int, TYPE*, int*, int*);

void ATL_CBC_Malloc_Lock(ATL_CBC_t*, ATL_INT);
void ATL_CBC_Malloc_UnLock(ATL_CBC_t*, ATL_INT);
/*
 * For updating the counter, a regular counter will casue overflow after a
 * few billion syncs. So, we will use boolean and toggle for each POST.
 */
/* void ATL_POST(enum ATL_SYNC_SCOPE, ATL_CBC_t*, ATL_INT); */
#define ATL_POST(scope, cbc, rankG)   (cbc->tdata[rankG]).Next[scope] ^= 1
#define ATL_ROW_RANK(cbc, rankG) ((cbc->tdata[rankG]).rankR)
#define ATL_COL_RANK(cbc, rankG) ((cbc->tdata[rankG]).rankC)
#define ATL_GRID_RANK(cbc, rankR, rankC) (cbc->Q*rankR + rankC)
#define ATL_NEXT(cbc, rankG, scope) ((cbc->tdata[rankG]).Next[scope])
#define ATL_TDATA(rankG) (cbc->tdata[rankG])

#endif
