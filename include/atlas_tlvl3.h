#ifndef atlas_tlvl3_H
   #define atlas_tlvl3_H
   #define DMM_H 1
   #define SMM_H 1
   #define CMM_H 1
   #define ZMM_H 1

#include "atlas_threads.h"
#ifdef TYPE
   #include "atlas_lvl3.h"
   #include "atlas_amm.h"
#endif
#ifndef ATL_XOVER_L3
   #ifdef TREAL
      #define ATL_XOVER_L3 2   /* number of NBxNB blocks */
   #else
      #define ATL_XOVER_L3 1
   #endif
#endif

#ifndef ATL_TGEMM_XOVER
   #define ATL_TGEMM_XOVER ATL_XOVER_L3
#endif
#ifndef ATL_TGEMM_ADDP
   #define ATL_TGEMM_ADDP 1
#endif
/*
 * Number of blocks per proc for GEMM to divide M only
 */
#ifndef ATL_TMMMINMBLKS
   #define ATL_TMMMINMBLKS 4
#endif
#ifndef ATL_TGEMM_THRESH_MF
   #define ATL_TGEMM_THRESH_MF \
      ((((2.0*(ATL_TGEMM_XOVER))*ATL_AMM_LLCMU)*ATL_AMM_LLCMU)*ATL_AMM_LLCMU)
#endif
/*
 * This is the minimal number of flops each thread requires once THRESH
 * is exceeded
 */
#ifndef ATL_TGEMM_PERTHR_MF
   #define ATL_TGEMM_PERTHR_MF \
      ((((2.0*ATL_TGEMM_ADDP)*ATL_AMM_LLCMU)*ATL_AMM_LLCMU)*ATL_AMM_LLCMU)
#endif
/*
 * For debugging, can define ATL_SERIAL_COMBINE, and then it any required
 * workspaces of C will be allocated before beginning parallel operations,
 * and all required combined will happen after parallel operations are
 * done.
 */
// #define ATL_SERIAL_COMBINE
#ifdef ATL_SERIAL_COMBINE
typedef struct ATL_CombNode ATL_combnode_t;
struct ATL_CombNode
{
   ATL_INT M, N, ldw, ldd;
   void *W, *D;                 /* Work and Destination */
   ATL_combnode_t *next;
};
#endif
/*
 * The array Cinfp holds C partitioning information.  This array holds a
 * list of pointers to nodes whose data I have not been able to combine
 * with my native C partition.  The first nCw entries contain the pointers
 * to the MMNODE of allocated C workspaces that I have not been able to
 * combine.  If my node has C in workspace, I am the first entry in this array.
 * Sometimes, a child thread has been combined with me that owned a piece of
 * the original C.  These values do not need to be combined (they were written
 * to the original C), but we need to combine the range of "owned" workspaces
 * so that we know when it is legal for a parent node to add into the space.
 * The final nCp entries of Cinfp entries of Cinfp hold these original pieces
 * that need to be combined to create larger owned partitions (starting from
 * the end of the array).  If the C ptr is NULL, that means that entry has
 * been subsumed into a new entry.
 */
typedef struct ATL_TMMNode ATL_TMMNODE_t;
struct ATL_TMMNode
{
   ATL_TMMNODE_t *Cinfp[ATL_NTHREADS];
   void (*gemmK)(ATL_CINT, ATL_CINT, ATL_CINT, const void*, const void *,
                 ATL_CINT, const void*, ATL_CINT, const void*, void*, ATL_CINT);
   const void *A, *B;
   void *C, *Cw;
   void *alpha, *beta;
   void *zero, *one;
   ATL_INT ldcw, M, N, K, lda, ldb, ldc;
   int mb, nb, kb;
   int eltsz, eltsh; /* element size, and shift (eg. log_2(eltsz)) */
   int rank;         /* the rank of my thread ([0,P-1]) */
   int nCw;          /* # of workspace entries in 1st nCw elts of Cinfp array */
   int nCp;          /* # of orig. C pieces last nCp elts of Cinfp */
   int ownC;         /* do I own my piece of C, or only wrkspace? */
};
/*
 * This data structure used for dynamically scheduled rank-K update
 * It is needed only by routines that are typed, and thus define TYPE
 */
#ifdef TYPE
typedef struct
{
   void *aNcnt;           /* count on col-panels of C */
   void *aMcnt;           /* count row-panels of A */
   void **aMcnts;         /* P-len array of counts on row-blks of C */
   void **Mlocks;         /* mutexes protecting init of aMcnts */
   int *Js;               /* current C col for each node */
   int Sync0;             /* 0: no sync at end; else thr 0 waits til all done */
   volatile int *chkin;   /* ATL_NTHREAD-len checkin array */
   TYPE **Bws;            /* preallocated thread copy areas */
   TYPE *Aw;              /* workspace for common A */
   const TYPE *A, *B;     /* original input matrices */
   TYPE *C;               /* original output matrix */
   #ifdef TREAL
      TYPE alpha;
      TYPE beta;
   #else
      const TYPE *alpha;
      const TYPE *beta;
   #endif
   ATL_INT nKb, kr, kr8;
   ATL_INT nMb, mr, nNb, nr;
   ATL_INT M, N, K, lda, ldb, ldc;
   enum ATLAS_TRANS TA, TB;
} ATL_TGEMM_RKK_t;
#endif

/*
 * This data structure is used when we split K for SYRK
 */
typedef struct ATL_SyrkK ATL_TSYRK_K_t;
struct ATL_SyrkK
{
   ATL_TSYRK_K_t *Cinfp[ATL_NTHREADS];
   void (*gemmT)(const enum ATLAS_TRANS, const enum ATLAS_TRANS,
                  ATL_CINT, ATL_CINT, ATL_CINT, const void *,
                  const void *, ATL_CINT, const void *, ATL_CINT,
                  const void *, void *, ATL_CINT);
   void (*tvsyrk)(const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CINT,
                  ATL_CINT, const void*, const void*, ATL_CINT, const void*,
                  void*, ATL_CINT);
   const void *A;
   void *C, *Cw;
   void *DoComb;
   ATL_LAUNCHSTRUCT_t *lp;
   const void *alpha, *beta;
   const void *zero, *one;
   ATL_INT ldcw, N, K, nb, lda, ldc;
   int eltsh, rank, nCw;
   enum ATLAS_UPLO Uplo;
   enum ATLAS_TRANS Trans, TB;
};
#ifdef TYPE
   #define VINT volatile int
/*
 * This structure used when K <= MAXKB, and M and N are large
 */
typedef struct ATL_tamm_rkK ATL_tamm_rkK_t;
struct ATL_tamm_rkK
{
   ammkern_t amm_b0;
   cm2am_t a2blk;       /* block copy for A */
   cm2am_t b2blk;       /* block copy for B, applies alpha */
   ablk2cmat_t blk2c;   /* copy that applies beta  */
   const TYPE *A;       /* input A matrix */
   const TYPE *B;       /* input B matrix */
   TYPE *C;             /* output matrix */
   TYPE *w;             /* nthr wsz-len thread-local workspaces */
   const TYPE *alpha;   /* ptr to alpha */
   const TYPE *beta;    /* ptr to beta */
   void *B1cpyAsgCtr;   /* 1: 1st B blk not assgnd, caller copies,0: cpy done */
   void *B1cpyDonCtr;   /* 1: 1st B blk not copied yet, 0: cpy done */
   void *AcpyCtr;       /* when 0, all of A has been copied */
   void **MbCtr;        /* nnblk-len array of Mblk ctrs */
   void *BAssgBV;       /* 1 means being copied, 0: not assigned */
   void *BDoneBV;       /* 1 means already copied, 0: not yet copied */
   void *cpBmut;        /* protects B[assg,done]BV */
   size_t wsz;          /* size of local workspace */
   int BCPYDONE;        /* if 1, all of B has been copied */
   int ACPYDONE;        /* if 1, all of A has been copied */
   int bsz;             /* size of common workspace for B */
   int TA;              /* 0: noTrans; 1: Trans */
   int TB;              /* 0: noTrans; 1: Trans */
   int N;               /* # cols of C; N <= MAXNB */
   int K;               /* common dim A&B; K <= MAXKB */
   int nmblks;          /* # of M blocks, including any partial block */
   int mr;              /* M%mu */
   int nbm;             /* # of M blocks 1 MbCtr gives out (1st col always 1) */
   int nmu;             /* CEIL(mb/mu) */
   int nmuL;            /* # of mus in final block */
   int mb;              /* block size for all blocks but last */
   int mbL;             /* block size for last block */
   int nnu;             /* CEIL(N/nu) */
   int KB0;             /* if kmajor, it is CEIL(K/ku)*ku, else K */
   int lda;             /* leading dim of A */
   int ldb;             /* leading dim of B */
   int ldc;             /* leading dim of C */
};
/*
 * This structure used when M <= MAXMB && N <= MAXNB (K large)
 */
typedef struct ATL_tamm_tMN ATL_tamm_tMN_t;
struct ATL_tamm_tMN
{
   ammkern_t amm_b0, amm_b1, ammK_b0, ammK_b1;
   cm2am_t a2blk;       /* block copy for A */
   cm2am_t b2blk;       /* block copy for B */
   const TYPE *A;       /* input A matrix */
   const TYPE *B;       /* input B matrix */
   TYPE *C;             /* output matrix */
   TYPE *w;             /* ATL_NTHREADS wsz-len thread-local workspaces */
   void *KbCtr;         /* which K block */
   int TA;              /* 0: noTrans; 1: Trans */
   int TB;              /* 0: noTrans; 1: Trans */
   int nC;              /* nmu*nnu*mu*nu */
   int M;               /* # of rows of C; M <= MAXMB */
   int N;               /* # cols of C; N <= MAXNB */
   int nmu;             /* CEIL(M/mu) */
   int nnu;             /* CEIL(N/nu) */
   int kb;              /* K blocking factor for all but first block */
   int kb0;             /* K blocking factor for first block */
   int KB0;             /* if k-vectorized, it is CEIL(kb0/ku)*ku, else  */
   int szA, szB;        /* mb*kb, nb*kb */
   size_t szW;          /* size of local workspace */
   size_t K;            /* common dim A&B; K is large & is parallelized */
   size_t nkblks;       /* # of K blocks, including any partial block */
   size_t lda;          /* leading dim of A */
   size_t ldb;          /* leading dim of B */
   size_t ldc;          /* leading dim of C */
};
/*
 * This structure used when N <= MAXNB && K <= MAXKB
 */
typedef struct ATL_tamm_tNK ATL_tamm_tNK_t;
struct ATL_tamm_tNK
{
   ammkern_t amm_b0;
   cm2am_t a2blk;       /* block copy for A */
   cm2am_t b2blk;       /* block copy for B, applies alpha */
   ablk2cmat_t blk2c;   /* copy for beta=1, applies beta  */
   const TYPE *A;       /* input A matrix */
   const TYPE *B;       /* input A matrix */
   TYPE *C;             /* output matrix */
   TYPE *w;             /* ATL_NTHREADS wsz-len thread-local workspaces */
   const TYPE *alpha;   /* ptr to alpha */
   const TYPE *beta;    /* ptr to beta */
   void *MbCtr;         /* which M block */
   void *BassgCtr;      /* 1 means must be copied, 0 assigned */
   void *BdoneCtr;      /* 0 means must be copied, 1 not ready */
   size_t wsz;          /* size of local workspace */
   int bsz;             /* size of common workspace for B */
   int TA;              /* 0: noTrans; 1: Trans */
   int TB;              /* 0: noTrans; 1: Trans */
   int N;               /* # cols of C; N <= MAXNB */
   int K;               /* common dim A&B; K <= MAXKB */
   int nmblks;          /* # of M blocks, including any partial block */
   int mr;              /* M%mu */
   int nmu;             /* CEIL(mb/mu) */
   int nmuL;            /* # of mus in final block */
   int mb;              /* block size for all blocks but last */
   int mbL;             /* block size for last block */
   int nnu;             /* CEIL(N/nu) */
   int KB0;             /* if kmajor, it is CEIL(K/ku)*ku, else K */
   int lda;             /* leading dim of A */
   int ldb;             /* leading dim of B */
   int ldc;             /* leading dim of C */
};
/*
 * This structure used when N <= MIN(MAXNB,MAXMB), so we deal out only
 * K blocks using the global counter KbCtr.
 */
typedef struct ATL_tsyrk_ammK ATL_tsyrk_ammK_t;
struct ATL_tsyrk_ammK
{
   ammkern_t amm_b0, amm_b1, ammK_b0, ammK_b1;
   cm2am_t a2blk;        /* no-transpose copy */
   cm2am_t b2blk;        /* transpose copy */
   ablk2cmat_t blk2c_b0;/* copy for beta=0 */
   ablk2cmat_t blk2c_b1;/* copy for beta=1 */
   const TYPE *A;       /* input matrix */
   TYPE *C;             /* output matrix */
   TYPE *w;             /* 2*ATL_NTHREADS thread-local mb*nb workspaces */
   const TYPE *alpha;   /* ptr to alpha */
   const TYPE *beta;    /* ptr to beta */
   void *KbCtr;         /* which k block */
   void *Cmut;          /* mutex lock for block of C */
   size_t wsz;          /* size of local workspace */
   VINT BETA_APPLIED;
   int LOWER;           /* set true if lower triangular C */
   int TA;              /* 0: noTrans; 1: Trans */
   int nkblks;          /* # of k blocks, including any partial block */
   int N;               /* total size of C, known to be <= MAXNB */
   int mb;              /* ((N+mu-1)/mu)*mu */
   int nb;              /* ((N+nu-1)/nu)*nu */
   int mbnb;            /* mb * nb */
   int nmu;             /* CEIL(N/mu) */
   int nnu;             /* CEIL(N/nu) */
   int kb;              /* blocking used for K */
   int kb0;             /* K%kb, if 0, kb */
   int KB0;             /* if kmajor, it is CEIL(kb0/ku)*ku, else kb0 */
   int lda;             /* leading dim of A */
   int ldc;             /* leading dim of C */
};
/*
 * This structure used by dynamic access-major SYRK
 * In the first phase, we work only on diagonal blocks, while copying both
 * A & A'.  For diag work, we parallelize both N & K dims so that the copy
 * is done as quickly as possible.  Threads coming in first choose differing
 * diag blks; diagonal blocks are dealt out cheaply using the dCtr global
 * counter (which starts at nnblks == ndiag).
 * Once all diagonal blocks are dealt out, new threads will start using
 * the atomic ctr array KbegCtr array to share K work for each diagonal.
 * both KbegCtr & KdonCtr are nnblk-len arrays of atomic counters.  Each
 * counter starts at nkblks.  Once the block pointed to by KbegCtr is
 * completely copied, the copying array increments the KdonCtr.  Only one
 * core per diag will get KdonCtr == 0 after doing his copy, and this
 * core will seize cdmut mutex in order to set the appropriate bit in
 * cpydonBV, which is a nnblks-length bit vector.  If the kth bit is set,
 * that means the ith row of A & jth col of A' has been copied.
 * Once a thread gets KbegCtr for a particular diag of 0, it means there's
 * no more work for this block of C, and so it will seize the appropriate
 * Cdmuts mutex which protects each diagonal block of C, and write its
 * finished contribution out to C.  The first such thread to ever seize
 * the mutex will scope dbetaBV to find this diagonal block needs beta applied;
 * while later threads will use beta=1.
 * Eventually, all diagonal work is finished, and the first processor to
 * get 0 for all dCtr & KbegCtr requests will set NODWORK=1, so later
 * threads don't have to query all the counters to know they should proceed
 * to non-diagonal work.
 *
 * For non-diagonal work, we count the number of non-diagonal blocks of C,
 * which is initially stored in the ncblks variable, which is protected
 * the cwmut mutex, which also protects the cblkBV, which is a ncblk-len BV.
 * A unset bit means that particular non-diagonal block has not yet been
 * assigned to a thread, while a 1 means it has.  The mutex also protects
 * the cpydone variable, which is set to 1 when cpydonBV has all bits set.
 * So, threads wanting to do non-diagonal work will find the first unset
 * bit in cblkBV, and then translate that to a (i,j) C block coordinate.
 * If the ith & jth bits are both set in cpydonBV (or cpydone is set), then
 * they will take that block as their own to do (in this phase, each thread
 * gets an individual block of C to do) by setting the bit in cblkBV.
 * When all bits are set in cblkBV, then all work has been dealt out, and
 * threads will exit once they say there is no more work to do.
 * The master process joins all created threads, and can delete data structures
 * safely after all joins succeed.
 */
typedef struct ATL_tsyrk_ammN ATL_tsyrk_ammN_t;
struct ATL_tsyrk_ammN
{
   ammkern_t amm_b0, amm_b1, ammK;
   cm2am_t a2blk;    /* no-transpose copy */
   cm2am_t b2blk;    /* transpose copy */
   ablk2cmat_t blk2d;/* copy for diagonal blocks */
   ablk2cmat_t blk2c;/* copy for non-diagonal blocks */
   const TYPE *A;    /* input matrix */
   TYPE *C;          /* output matrix */
   TYPE *wC;         /* ATL_NTHREADS thread-local nb*nb workspaces */
   TYPE *wA, *wAt;   /* workspaces for storing A & A' */
   const TYPE *alpha;/* ptr to alpha */
   const TYPE *beta; /* ptr to beta */
   int *cpydonBV;    /* nnblks BV set means A & A' copy is done */
   int *cblkBV;      /* ncblks BV for dealing out cblks */
   int *dbetaBV;     /* nnblks BV: unset means beta not yet applied */
   int *cbetaBV;     /* ncblks BV: unset means beta not yet applied */
   void *cdmut;      /* mutex protecting cpydonBV */
   void *cwmut;      /* mutex protecting non-diag work */
   void *dCtr;       /* AtomicCtr for dealing out diagonal blocks */
   void **KbegCtr;   /* K-counters for dealing out diagonal blocks */
   void **KdonCtr;   /* K-counters for completed diagonal blocks */
   void **Cdmuts;    /* mutex locks for each diagonal block of C */
   size_t panszA;    /* nkblks * blkszA */
   #if defined(ATL_PHI_SORTED) || defined(ATL_PHI_SORTED2)
      VINT *chkin;   /* ncores*ATL_Cachelen check-in array */
      int ncores;    /* ncores, on PHI it is nthr/4 */
      int ncntxts;   /* ncontexts to use per core */
   #endif
   #ifdef ATL_PHI_SORTED2
      int nnuC;
      int CoffC;
      int BoffC;
   #endif
   VINT cpydone;     /* set when all A/A' copying complete */
   int LOWER;        /* set true if lower triangular C */
   int TA;           /* 0: noTrans; 1: Trans */
   int ndiag;        /* # of diagonal blocks */
   int ncblks;       /* # of non-diagonal blocks left to be assigned */
   int nkblks;       /* # of k blocks, including any partial block */
   int NODWORK;      /* set to 1 by first thread to find all work started */
   int nb;           /* N blocking used for all but remainder block */
   int nbnb;         /* nb * nb */
   int nmu;          /* nb/mu */
   int nnu;          /* nb/nu */
   int nbf;          /* N%nb, if 0, nb */
   int nnuf;         /* CEIL(nbf/nu) */
   int nmuf;         /* CEIL(nbf/mu) */
   int kb;           /* blocking used for K */
   int kb0;          /* K%kb, if 0, kb */
   int KB0;          /* if kmajor, it is CEIL(kb0/ku)*ku, else kb0 */
   int blkszA;       /* nb*kb */
   int Mf;           /* ((nbf+mu-1)/mu)*mu */
   int Nf;           /* ((nbf+nu-1)/nu)*nu */
   int lda;          /* leading dim of A */
   int ldc;          /* leading dim of C */
};
/*
 * This data structure for doing access-major threaded gemm for moderately
 * sized problems where no dimension is <= to the blocking factor, and
 * the problem is not too large to prevent us from copying all of A & B
 * up front.  Copying up front allows us to compute the blocks of C
 * in any order.  For large problems, will have to recur (mainly on K) to
 * get workspace down to reasonable levels.  This case will not work well
 * if the number of C blocks is not quite a bit larger than the nthreads.
 */

typedef struct ATL_tgemm_ammG ATL_tgemm_ammG_t;
struct ATL_tgemm_ammG
{
   ammkern_t amm_b0, amm_b1, ammK_b0, ammK_b1;
   cm2am_t a2blk;    /* copy/transpose for A */
   cm2am_t b2blk;    /* copy/tranpose for B */
   ablk2cmat_t blk2c,/* access-major to column-major copy/scale for C */
     blk2c_b1;       /* access-major to column-major copy for C, BETA=1 */
   const TYPE *A;    /* input matrix */
   const TYPE *B;    /* input matrix */
   TYPE *C;          /* output matrix */
   TYPE *wC;         /* ATL_NTHREADS thread-local mb*nb workspaces */
   TYPE *wA, *wB;    /* workspaces for storing A & B */
   SCALAR beta;      /* beta */
   SCALAR alpA;      /* alpha to apply to A */
   SCALAR alpB;      /* alpha to apply to B */
   SCALAR alpC;      /* alpha to apply to C */
   int *cpyAdBV;     /* nmblks BV, set means K-panel of A has been copied */
   int *cpyBdBV;     /* nnblks BV, set means K-panel of A has been copied */
   int *cCblkBV;     /* nMNblks BV for dealing C blks while copying */
   int *cbetaBV;     /* nMNblks BV: unset means beta not yet applied */
   void *cpmut;      /* mutex protecting cpyA/BdBV & cCblkBV */
   void *cbetamut;   /* mutex protecting cbetaBV */
   void *ccCtr;      /* nMNblks ctr for dealing out (diagonal) blocks wt copy */
   void *cCtr;       /* ncblks ctr for dealing out blocks w/o copy */
   void **KbegCtr;   /* counters for dealing out K blocks for copying */
   void **KdonCtr;   /* K-counters for copied blocks */
   void **Cmuts;     /* MNblks mutex locks for copy-blocks of C */
   size_t panszA;    /* nkblks * blkszA */
   size_t panszB;    /* nkblks * blkszB */
   size_t lda;       /* leading dim of A */
   size_t ldb;       /* leading dim of B */
   size_t ldc;       /* leading dim of C */
   #ifdef ATL_PHI_SORTED
      VINT *chkin;   /* ncores*ATL_Cachelen check-in array */
      int ncores;    /* ncores, on PHI it is nthr/4 */
      int ncntxts;   /* ncontexts to use per core */
   #endif
   VINT cpyAdone;    /* set when all of A has been copied */
   VINT cpyBdone;    /* set when all of B has been copied */
   VINT NOCPWORK;    /* set to 1 by 1st thread to find all copy work started */
   int TA;           /* 0: noTrans; 1: Trans */
   int TB;           /* 0: noTrans; 1: Trans */
   int nCblks;       /* nmblks * nnblks */
   int nMNblks;      /* MAX(nmblks,nnblks) */
   int nmblks;       /* CEIL(M/mb) */
   int nnblks;       /* CEIL(N/nb) */
   int nkblks;       /* CEIL(K/kb) */
   int mb;           /* M blocking used for all but remainder block */
   int nb;           /* N blocking used for all but remainder block */
   int nmu;          /* mb/mu */
   int nnu;          /* nb/nu */
   int mbf;          /* M%mb, if 0, mb */
   int nbf;          /* N%nb, if 0, nb */
   int nmuf;         /* CEIL(mbf/mu) */
   int nnuf;         /* CEIL(nbf/nu) */
   int kb;           /* blocking used for K */
   int kb0;          /* K%kb, if 0, kb */
   int KB0;          /* if kmajor, it is CEIL(kb0/ku)*ku, else kb0 */
   int blkszA;       /* mb*kb */
   int blkszB;       /* nb*kb */
   int blkszC;       /* mb*nb */
};
   #undef VINT
#endif
/*
 * This data structure used when we divide N only, and NTHREAD is a power of 2
 */
typedef struct
{
   void (*gemmK)(ATL_CINT, ATL_CINT, ATL_CINT, const void*, const void *,
                 ATL_CINT, const void*, ATL_CINT, const void*, void*, ATL_CINT);
   void (*tvsyrk)(const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CINT,
                  ATL_CINT, const void*, const void*, ATL_CINT, const void*,
                  void*, ATL_CINT);
   void *T;             /* Triangular matrix to do SYRK into*/
   void *C;             /* rect matrix to do GEMM into */
   const void *A0;      /* input matrix for syrk, */
   const void *A;       /* 1st input matrix for GEMM */
   const void *B;       /* 2nd input matrix for GEMM */
   const void *alpha, *beta;
   ATL_INT M;           /* size of SYRK and 1st dim of GEMM */
   ATL_INT N;           /* size of 2nd dim of N */
   ATL_INT K;           /* K of original problem */
   ATL_INT lda, ldc;
   int nb, eltsh;       /* shift to do equivalant of *sizeof */
   enum ATLAS_UPLO Uplo;
   enum ATLAS_TRANS TA, TB;
} ATL_TSYRK_M_t;


typedef struct
{
   const void *A, *alpha;
   void *B;
   ATL_INT M, N, lda, ldb;
   enum ATLAS_SIDE side;
   enum ATLAS_UPLO uplo;
   enum ATLAS_TRANS TA;
   enum ATLAS_DIAG  diag;
} ATL_TTRSM_t;

typedef struct
{
   const void *A, *B, *alpha, *beta;
   void *C;
   ATL_INT M, N, lda, ldb, ldc, nb;
   enum ATLAS_SIDE side;
   enum ATLAS_UPLO uplo;
} ATL_TSYMM_t;
typedef struct
{
   const void *alpha, *alpha2, *beta, *one, *zero;
   void (*tvgemm)(const enum ATLAS_TRANS, const enum ATLAS_TRANS,
                  ATL_CINT, ATL_CINT, ATL_CINT, const void *,
                  const void *, ATL_CINT, const void *, ATL_CINT,
                  const void *, void *, ATL_CINT);
   void (*tvApAt)(const enum ATLAS_UPLO, ATL_CINT, const void *, ATL_CINT,
                  const void *, void *, ATL_CINT);

   ATL_INT K, lda, ldb, ldc;
   int nb, eltsh;
   enum ATLAS_UPLO Uplo;
   enum ATLAS_TRANS trans, TA, TB,  /* trans for syr2k, TA,TB are for GEMM */
                    TA2, TB2;       /* transpose of TA,TB */
} ATL_SYR2K_t;

/*
 * For triangular matrices, diagonal blocks are handled specially, but we
 * get dense square blocks above/below the diagonal.  We consider Upper
 * triangular the transpose of Lower, allowing us to only handle Lower.
 * Our AtomicCtr routines are 1-D counters, not 2-D, so we linearize the
 * blocks beneath the diagonal by counting them column-wise.  So, a 4x4
 * matrix of blocks would look like:
 *    |X X X X|
 *    |1 X X X|
 *    |2 4 X X|
 *    |3 5 7 X|
 */
/*
 * Translates (i,j) coordinate of lower triangular matrix to number between
 * [0,n), where n = number of non-diagonal block.  nm_ is the number of
 * diagonal blocks in the matrix (4 in example above).
 */
#define Mcoord2tblk(nm_, i_, j_) \
   ((nm_)*(j_) - (((1+(j_))*(j_))>>1) + (i_)-(j_)-1)

/*
 * In a lower matrix of with NB_ diagonal blocks, translate the linearized
 * block number B_ back into rectangular (I,J) coordinates.  The difficulty
 * is finding J.  Would like to do it with an equation, like we do when
 * converting from coordinates to block number.  Tony Castaldo came up with
 *    J = (int)((nm-0.5-0.5*sqrt(4*nm*nm-4*nm+1-8*b)))
 * but sqrt is a function call which does a Newtonian iteration on floats
 * (therefore, has a relatively slow loop).  Finding the J column is indeed
 * the hard part, and in theory we can use binary search to find in
 * O(log_2(NB_)) time.  However, this algorithm requires multiplication
 * inside the loop, and so it is never competitive for the range we are
 * interested in (NB_ usually < 10, and almost never larger than 2000).
 * So instead do a linear search to find j, but optimize by first
 * constraining j to a 128-column region, then a 8-column region, and
 * then find the actual column.  So worst-case loop counts are
 *  CEIL(NB_/128) + 16 + 8
 *
 * nblkC = # of blocks in a [128,8]-column Chunk
 */
#define Mtblk2coord(NB_, B_, I_, J_) \
{ \
   unsigned int j_=0; \
   unsigned int n_ = (NB_), b_=(B_), nblksC_; \
   KEEP_LOOKING128: /* find 128-col region where J is */ \
      nblksC_ = (n_<<7)-8256; \
      if (b_ < nblksC_ || n_ < 128) \
         goto FOUND128; \
      b_ -= nblksC_; \
      n_ -= 128; \
      j_ += 128; \
   goto KEEP_LOOKING128; \
   FOUND128: \
   KEEP_LOOKING8: /* find 8-col region where J is */ \
      nblksC_ = (n_<<3)-36; \
      if (b_ < nblksC_ || n_ < 8) \
         goto FOUND8; \
      b_ -= nblksC_; \
      n_ -= 8; \
      j_ += 8; \
   goto KEEP_LOOKING8; \
   FOUND8: \
      for (n_--; b_ >= n_; j_++) \
         b_ -= n_--; \
   (J_) = j_; \
   (I_) = j_ + b_ + 1; \
}

/*
 * =============================================================================
 * Function prototypes
 * =============================================================================
 */
void ATL_EnforceNonPwr2LO(ATL_TMMNODE_t *ptmms, const int P);
int Mjoin(PATL,threadMM)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                         size_t M, size_t N, size_t K);
void ATL_tvsyr2k_rec(ATL_SYR2K_t *syp, ATL_CINT Nblks, ATL_CINT nr,
                     const void *A, const void *B, void *C);
#ifdef TYPE
void Mjoin(PATL,tsyrk)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans, ATL_CINT N,
    ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
#ifdef TCPLX
void Mjoin(PATL,therk)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans, ATL_CINT N,
    ATL_CINT K, const TYPE alpha, const TYPE *A, ATL_CINT lda,
    const TYPE beta, TYPE *C, ATL_CINT ldc);
#endif
void Mjoin(PATL,tsymm)
   (const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc);
#ifdef TCPLX
void Mjoin(PATL,themm)
   (const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc);
#endif
void Mjoin(PATL,tsyr2k)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsyr2k)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc);
#ifdef TCPLX
void Mjoin(PATL,ther2k)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const TYPE beta, TYPE *C, ATL_CINT ldc);
#endif

void Mjoin(PATL,ttrsm)(const enum ATLAS_SIDE side, const enum ATLAS_UPLO uplo,
                       const enum ATLAS_TRANS TA, const enum ATLAS_DIAG diag,
                       ATL_CINT M, ATL_CINT N, const SCALAR alpha,
                       const TYPE *A, ATL_CINT lda, TYPE *B, ATL_CINT ldb);
void Mjoin(PATL,ttrmm)(const enum ATLAS_SIDE side, const enum ATLAS_UPLO uplo,
                       const enum ATLAS_TRANS TA, const enum ATLAS_DIAG diag,
                       ATL_CINT M, ATL_CINT N, const SCALAR alpha,
                       const TYPE *A, ATL_CINT lda, TYPE *B, ATL_CINT ldb);
void Mjoin(PATL,tgemm)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                       ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
                       const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
                       const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tvgemm)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                        ATL_CINT M, ATL_CINT N, ATL_CINT K, const void *alpha,
                        const void *A, ATL_CINT lda, const void *B,ATL_CINT ldb,
                        const void *beta, void *C, ATL_CINT ldc);
int Mjoin(PATL,ammm_REC)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc,
    int (*amm)(enum ATLAS_TRANS,enum ATLAS_TRANS, ATL_CINT, ATL_CINT, ATL_CINT,
               const SCALAR, const TYPE*, ATL_CINT,  const TYPE*, ATL_CINT,
               const SCALAR, TYPE*, ATL_CINT));
int Mjoin(PATL,tammm_tMN)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
int Mjoin(PATL,tammm_tNK)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
int Mjoin(PATL,tammm_G)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);

int Mjoin(PATL,tgemm_bigMN_Kp)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
int Mjoin(PATL,tgemm_rkK)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
int Mjoin(PATL,tgemm_rec)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
#ifdef TCPLX
void Mjoin(PATL,tgemmCC)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmCC)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
#endif  /* end ifdef TCPLX */
#ifdef TCPLX
void Mjoin(PATL,tgemmCN)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmCN)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
#endif  /* end ifdef TCPLX */
#ifdef TCPLX
void Mjoin(PATL,tgemmCT)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmCT)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
#endif  /* end ifdef TCPLX */
#ifdef TCPLX
void Mjoin(PATL,tgemmNC)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmNC)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
#endif  /* end ifdef TCPLX */
void Mjoin(PATL,tgemmNN)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmNN)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
void Mjoin(PATL,tgemmNT)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmNT)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
#ifdef TCPLX
void Mjoin(PATL,tgemmTC)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmTC)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
#endif  /* end ifdef TCPLX */
void Mjoin(PATL,tgemmTN)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmTN)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
void Mjoin(PATL,tgemmTT)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tsvgemmTT)
   (ATL_CINT M, ATL_CINT N, ATL_CINT K, const void* alpha,
    const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
    const void *beta, void *C, ATL_CINT ldc);
#endif  /* end ifdef TYPE */

#ifdef ATL_ESTNCTR
   #ifdef __GNUC__
   static inline int ATL_EstNctr    /* RETURNS: good # of global ctrs */
   #else
   static int ATL_EstNctr           /* RETURNS: good # of global ctrs */
   #endif
   (
      int N,   /* Total count you are going to use */
      int P    /* nthreads to be used */
   )
   {
      int ncnt, ncntP;
/*
 *    Want at least 32 elements per counter
 */
      ncnt = (N > 64) ? (N >> 5) : 1;
/*
 *    allow up to 8-way contention
 */
      if (P >= 16)
         ncntP = (P>>3);
      else
         ncntP = (P >= 4) ? (P>>1) : 1;
      return(Mmin(ncnt,ncntP));
   }
#endif

#endif
