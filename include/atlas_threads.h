#ifndef ATLAS_THREADS_H
   #define ATLAS_THREADS_H
#include "atlas_ttypes.h"
#include "atlas_taffinity.h"
/*
 * Unless told otherwise, use the less-intrusive conditional variable-based
 * thread pool, unless we're on the PHI, where that scale and slow OS make
 * that prohibitive.
 * NOTE: right now, always use sleep, as there may be a race condition in
 *       polling thread pool.
 */
#if !defined(ATL_POLLTPOOL) && !defined(ATL_SLEEPTPOOL)
   #define ATL_SLEEPTPOOL 1
#endif
/*
 * Need to fix this later by using a probe.  ATL_PHI_SORTED being defined
 * asserts that the first P/4 cores in the main list use context 0, the
 * next context 1, and so on.  This is the default sorting by ATLAS, and
 * we'll presently assume it true as long as the NTHR%4 == 0
 */
#if !defined(ATL_PHI_SORTED) && !defined(ATL_PHI_UNSORTED)
   #if defined(ATL_ARCH_XeonPHI) && ATL_NTHREADS%4 == 0
      #define ATL_PHI_SORTED 1
   #else
      #define ATL_PHI_UNSORTED 1
   #endif
#endif
#ifndef ATL_PTMAXMALLOC
   #ifndef ATL_PTMAXMALLOC_MB
      #if defined(ATL_ARCH_XeonPHI)
         #define ATL_PTMAXMALLOC_MB 32UL
      #else
         #define ATL_PTMAXMALLOC_MB 128UL
      #endif
   #endif
   #define ATL_PTMAXMALLOC (((ATL_PTMAXMALLOC_MB)<<20)<<ATL_NTHRPOW2)
   #ifdef ATL_MaxMalloc
      #if ATL_MaxMalloc*ATL_NTHREADS < ATL_PTMAXMALLOC
         #undef ATL_PTMAXMALLOC
         #define ATL_PTMAXMALLOC (ATL_MaxMalloc*ATL_NTHREADS)
      #elif ATL_MaxMalloc > ATL_PTMAXMALLOC
         #undef ATL_PTMAXMALLOC
         #define ATL_PTMAXMALLOC ATL_MaxMalloc
      #endif
   #endif
#endif
/*
 * If we don't have thread affinity, then the thread I'm waiting on may share
 * the core with me.  In this case, yield my time slice when polling
 */
#include "atlas_tsumm.h"
#include <limits.h>
#if defined(ATL_TAFFINITY) && ATL_TAFFINITY
   #define ATL_POLL
#else
   #define ATL_POLL ATL_thread_yield()
#endif

#if defined(ATL_OS_Win64) || defined(ATL_OS_WinNT)
   #ifdef ATL_USE64BITS
      #define ATL_WIN64THREADS 1
      #define ATL_WINTHREADS 1
   #else
      #define ATL_WIN32THREADS 1
      #define ATL_WINTHREADS 1
   #endif
#endif
#include "atlas_pthreads.h" /* gened file defs ATL_NTHREADS & ATL_NTHRPOW2 */
#ifdef ATL_WINTHREADS
   #include <windows.h>
   typedef struct
   {
      void *vp;      /* ptr to extra info */
      HANDLE thrH;   /* handle to thread */
      int rank;      /* my rank */
      int P;         /* # of processors in this call */
      int affID;     /* affinity id for this core */
   } ATL_thread_t;
   #ifndef CREATE_SUSPENDED
      #define CREATE_SUSPENDED 0x00000004
   #endif
   #ifndef WAIT_FAILED
      #define WAIT_FAILED (~0)
   #endif
#elif defined(ATL_OMP_THREADS)
   #include <omp.h>
   typedef struct
   {
      void *vp;      /* ptr to extra info */
      int rank;      /* my rank */
      int P;         /* # of processors in this call */
      int affID;     /* < 0: not set, affID=1-affID, else affID */
   } ATL_thread_t;
#else
   #define ATL_USE_POSIXTHREADS 1
   #include <pthread.h>
   typedef struct
   {
      pthread_t thrH;/* handle of thread */
      void *vp;      /* ptr to extra info */
      int rank;      /* my rank */
      int P;         /* # of processors in this call */
      int affID;     /* < 0: not set, affID=1-affID, else affID */
   } ATL_thread_t;
#endif

typedef struct ATL_LaunchStruct ATL_LAUNCHSTRUCT_t;
struct ATL_LaunchStruct
{
   ATL_thread_t *rank2thr;              /* index by rank to get thread handle */
   void *opstruct;
   int (*OpStructIsInit)(void*);        /* Query to see if struct is valid */
   void (*DoWork)(ATL_LAUNCHSTRUCT_t*, void*);
   void (*DoComb)(void*, const int, const int);  /* combine func */
   void *vp;                            /* misc. extra info ptr */
   volatile int *chkin;                 /* nthr-len checkin array */
   void **acounts;                      /* var-len array of atomic counters */
};

/*
 * Constants for use in chkin array
 */
#define ATL_CHK_DONE_OP   -2048
#define ATL_CHK_DONE_COMB -2049

/*
 * The following info is all for when a thread pool is being used, rather
 * than launch & join paradigm.
 */
/*
 * Thread pool flag macros
 */
#define ATL_TPF_POLL(p_) ((p_)->pflag & 1)
#define ATL_TPF_SET_POLL(p_) ((p_)->pflag |= 1)
#define ATL_TPF_UNSET_POLL(p_) ((p_)->pflag &= ~1)
#define ATL_TPF_ZEROWAKES(p_) ((p_)->pflag & 2)
#define ATL_TPF_SET_ZEROWAKES(p_) ((p_)->pflag |= 2)
#define ATL_TPF_UNSET_ZEROWAKES(p_) ((p_)->pflag &= ~2)
#define ATL_TPF_DYNCOMB(p_) ((p_)->pflag & 4)
#define ATL_TPF_SET_DYNCOMB(p_) ((p_)->pflag |= 4)
#define ATL_TPF_UNSET_DYNCOMB(p_) ((p_)->pflag &= ~4)
#define ATL_TPF_DIE(p_) ((p_)->pflag & 8)
#define ATL_TPF_SET_DIE(p_) ((p_)->pflag |= 8)
#define ATL_TPF_UNSET_DIE(p_) ((p_)->pflag &= ~8)
#define ATL_TPF_MICSORTED(p_) ((p_)->pflag & 16)
#define ATL_TPF_SET_MICSORTED(p_) ((p_)->pflag |= 16)
#define ATL_TPF_UNSET_MICSORTED(p_) ((p_)->pflag &= ~16)

#define VUINT volatile unsigned int
/*
 * Presently, always-polling only threadpool option designed to work on Windows
 */
#if defined(ATL_WINTHREADS)
   #define ATL_TP_FULL_POLL
#elif defined(ATL_SLEEPTPOOL) && ATL_SLEEPTPOOL
   #ifdef ATL_TP_FULLPOLL
      #undef ATL_TP_FULLPOLL
   #endif
#elif defined(ATL_TP_FULLPOLL) && !ATL_TP_FULLPOLL
   #undef ATL_TP_FULLPOLL
#elif defined(ATL_HAS_AFFINITY) && ATL_HAS_AFFINITY && !defined(ATL_OMP_THREADS)
   #define ATL_TP_FULLPOLL 1
#endif
/*
 * Function pointer taking the pool as an argument that does the work.
 */
typedef void (*ATL_tpjfunc_t)(void *pp, int rank, int vrank);
/*
 * A thread pool takes can combine results at end of run if needed.
 * It takes everything the jobfunc does, + the vranks to combine
 */
typedef void (*ATL_tpjcomb_t)
   (void *pp, int rank, int vrank, int hisvrank);
/*
 * Definition of an ATLAS thread pool
 */
typedef struct ATL_ThreadPool ATL_tpool_t;
#ifdef ATL_TP_FULLPOLL
struct ATL_ThreadPool
{
   VUINT nthr;            /* # of threads total in this thread pool */
   VUINT pflag;           /* bitvector of options */
   VUINT nworkers;        /* # of thr supposed to wake up and work */
   #ifdef ATL_PHI_SORTED
      VUINT ncores, cntxts, nwcores;
   #endif
   volatile short *chkin; /* nthr-len chkin array wt meaning:
      /* 0: keep polling, 1: start, 2: have started, 3: finished 4: exited */
      /* 0&1 written by master, 2-4 written by worker */
   volatile int *icomm;   /* nthr-long integer communication array */
   ATL_thread_t *threads; /* array of thread ptrs */
/*
 * Variables below set only when doing a combine, and protected by tpmut
 * (which need only be set when doing a combine)
 */
   void *combmut;         /* mutex for doing optional combine op */
   int *combReadyBV;      /* 1: thr of that rank is ready to do combine */
   int *combDoneBV;       /* 1: thr of that rank's data already combined */
   void *tpmut;           /* mutex protecting above pool info */
/*
 * variables below here manipulated only when threads known to be polling
 */
   ATL_tpjfunc_t jobf;    /* ptr to job function to execute */
   ATL_tpjcomb_t combf;   /* NULL: don't combine, else combine func */
   void *PD;              /* problem def to give to jobf & combf */
   void *extra;           /* extra info for jobf & combf */
};
#elif defined(ATL_OMP_THREADS)
struct ATL_ThreadPool
{
   VUINT nthr;            /* # of threads total in this thread pool */
   VUINT pflag;           /* bitvector of options */
   VUINT nworkers;        /* # of thr supposed to wake up and work */
   #ifdef ATL_PHI_SORTED
      VUINT ncores, cntxts, nwcores;
   #endif
   volatile int *icomm;   /* nthr-long integer communication array */
   ATL_thread_t *threads; /* array of thread ptrs */
/*
 * Variables below set only when doing a combine, and protected by tpmut
 * (which need only be set when doing a combine)
 */
   void *combmut;         /* mutex for doing optional combine op */
   int *combReadyBV;      /* 1: thr of that rank is ready to do combine */
   int *combDoneBV;       /* 1: thr of that rank's data already combined */
   void *tpmut;           /* mutex protecting above pool info */
/*
 * variables below here manipulated only when threads known to be polling
 */
   ATL_tpjfunc_t jobf;    /* ptr to job function to execute */
   ATL_tpjcomb_t combf;   /* NULL: don't combine, else combine func */
   void *PD;              /* problem def to give to jobf & combf */
   void *extra;           /* extra info for jobf & combf */
};
#else
struct ATL_ThreadPool
{
   VUINT jobID;           /* count of jobs, wraps at top of uint range */
   VUINT WORKDONE;        /* zeroed to start job, set by last worker done */
   VUINT NOWORK;          /* optionally set when all work dealt out */
   VUINT nthr;            /* # of threads total in this thread pool */
   VUINT nsleep;          /* # of threads that have gone to sleep at start */
   VUINT nworkers;        /* # of thr supposed to wake up and work */
   VUINT wcnt;            /* count incremented as workers wake up */
   VUINT nwdone;          /* # of workers that have completed the task */
   VUINT pflag;           /* bitvector of options */
   void *wcond;           /* cond var for work pool sleep/wake */
   #ifdef ATL_PHI_SORTED
      void *wcond2;       /* cond vars 2nd context sleeps on */
      void *wcond3;       /* cond vars 3rd context sleeps on */
      void *wcond4;       /* cond vars 4th context sleeps on */
   #elif defined(ATL_PHI_SORTED2)
      int ncntxts;
      void **wconds;
   #endif
   void *mcond;           /* cond for master process sleep/wake */
   ATL_thread_t *threads; /* array of thread ptrs */
   volatile int *icomm;   /* nthr-long integer communication array */
   void *combmut;         /* mutex for doing optional combine op */
   int *combReadyBV;      /* 1: thr of that rank is ready to do combine */
   int *combDoneBV;       /* 1: thr of that rank's data already combined */
   void *tpmut;           /* mutex protecting above pool info */
/*
 * variables below here manipulated only when threads known asleep
 */
   ATL_tpjfunc_t jobf;    /* ptr to job function to execute */
   ATL_tpjcomb_t combf;   /* NULL: don't combine, else combine func */
   void *PD;              /* problem def to give to jobf & combf */
   void *extra;           /* extra info for jobf & combf */
};
#endif

/*
 * Declare the beautiful global variables used by thread pool
 */
#ifdef ATL_TP_DECL
   double ATL_POLLTIME=0.01;  /* poll for 10 millisecond */
   ATL_tpool_t *ATL_TP_PTR=NULL, *ATL_TP1_PTR=NULL;
#else
   extern double ATL_POLLTIME;
   extern ATL_tpool_t *ATL_TP_PTR, *ATL_TP1_PTR;
#endif
void *ATL_threadpool(void *tp);
void *ATL_threadpool_launch(void *tp);
void ATL_InitThreadPoolStartup(int P, void *pd, void *extra);
void ATL_goParallel (const unsigned int P, void *DoWork, void *DoComb,
                     void *PD, void *extra);
ATL_tpool_t *ATL_NewThreadPool(const int P, int ICOM, void *vp);
void ATL_FreeThreadPool(ATL_tpool_t *pp);
void ATL_oldjobwrap(void *vpp, int rank, int vrank);
void ATL_oldcombwrap(void *vpp, int rank, int vrank, int vhisrank);
int ATL_tpool_dojob(ATL_tpool_t *pp, const int rank, const int CallFrWrk);
void *ATL_threadpool(void *vp);
#undef VUINT
int ATL_setmyaffinity(const int); /* sets affinity of already-running thread */
/* Sets up ATL_LAUNCHSTRUCT_t var and ATL_thread_t array & starts threads*/
void ATL_thread_launch(void *opstruct, int opstructstride, void *OpStructIsInit,
                       void *DoWork, void *CombineOpStructs);
void ATL_goparallel(const unsigned int P, void *DoWork, void *opstruct, void*);
void ATL_goparallel_prank(const unsigned int P, void *DoWork, void *opstruct,
                          void *DoComb);
/*  Starts a thread running, and sets its affinity to proc if possible */
int ATL_thread_start(ATL_thread_t *thr, int proc, int JOINABLE,
                     void *(*rout)(void*), void*);
int ATL_thread_join(ATL_thread_t*); /* waits on completion of thread */
void ATL_thread_exit(void*);        /* exits currently executing thread */
void *ATL_log2tlaunch(void *vp);    /* min spanning tree launch */
void *ATL_lin0tlaunch(void *vp);    /* 0 linear launches all threads */
void *ATL_dyntlaunch(void *vp);     /* launch done as workpool */
/*
 * Atomic count functions; may be less overhead than mutex on some systems
 */
void *ATL_SetAtomicCount(int cnt);   /* allocates acnt, sets acnt=cnt */
int   ATL_ResetAtomicCount(void *vp, int cnt);  /* reset vp to cnt */
int   ATL_DecAtomicCount(void *vp);  /* returns acnt-- (not --acnt!) */
int   ATL_GetAtomicCount(void *vp);  /* returns acnt */
void  ATL_FreeAtomicCount(void *vp); /* free acount resources */
/*
 * Global count functions, built out of P local counters for good scaling
 */
void *ATL_SetGlobalAtomicCount(int P, int cnt, int percLoc);
void  ATL_ResetGlobalAtomicCount(void *vp, int cnt, int percLoc);
int   ATL_DecGlobalAtomicCount(void *vp, int rank);
int   ATL_GetGlobalAtomicCount(void *vp, int rank);
void  ATL_FreeGlobalAtomicCount(void *vp);
/*
 * Countdown funcs: same as above Glob, but guarantee 1 is last non-zero # ret
 */
void *ATL_SetGlobalAtomicCountDown(int P, int cnt);
int ATL_DecGlobalAtomicCountDown(void *vp, int rank);
void ATL_FreeGlobalAtomicCountDown(void *vp);

/*
 * If using pthreads, just wrapper around pthread mutex funcs, else written
 * in terms of the AtomicCount funcs with init value set to 1
 */
void *ATL_mutex_init(void);       /* returns mutex pointer */
void ATL_mutex_free(void *vp);    /* frees mutex vp */
void ATL_mutex_lock(void *vp);
void ATL_mutex_unlock(void *vp);
int  ATL_mutex_trylock(void *vp); /* opp pthreads: 0 means lock NOT acquired */
void ATL_thread_yield(void);      /* gives up time slice */
/*
 * Condition variables only used for thread pool, not implemented yet on
 * anything but pthreads (Windows & OpenMP missing)
 */
void *ATL_cond_init(void);
void ATL_cond_free(void *vp);
void ATL_cond_signal(void *cond);
void ATL_cond_bcast(void *cond);
void ATL_cond_wait(void *cond, void *mut);

#define MindxT(A_,i_) ((void*)( ((char*)(A_)) + ((size_t)i_) ))
#define ATL_tlaunch ATL_log2tlaunch   /* may want linear launch later */
void ATL_tDistMemTouch(size_t N, void *vp);

#endif   /* end of #ifdef protecting include file from redundant inclusion */

