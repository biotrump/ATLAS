/*
 * This file defines the types for the new threaded routines, without using
 * any specific thread info so that they can be safely included prior to
 * building & tuning the threading system
 */
#ifdef TYPE
#ifndef ATLAS_TTYPES_H
   #define ATLAS_TTYPES_H
#include "atlas_amm.h"
typedef volatile int VINT;
/*
 * This structure is for a TRSM that uses ATL_trsm as its compute kernel,
 * and only divides N (RHS) matrix.  Works best for tiny TRiangular matrix.
 */
typedef struct ATL_ttrsm_tTR_t ATL_ttrsm_tTR_t;
struct ATL_ttrsm_tTR_t
{
   void *rhsBlkCtr;        /* deals out RHS, [1-nrblks] */
   size_t lda, ldb;        /* leading dims of A & B */
   const TYPE *A;          /* triangular matrix */
   TYPE *B;                /* MxN input/output matrix of right-hand sides */
   #ifdef TREAL
      TYPE alpha;          /* scalar for B */
   #else
      const TYPE *alpha;   /* scalar for B */
   #endif
   enum ATLAS_SIDE side;   /* Whether B on Left or Right of A */
   enum ATLAS_TRANS TA;    /* Whether A is transposed or not */
   enum ATLAS_DIAG uplo;   /* Whether triangle stored Upper of Lower */
   enum ATLAS_DIAG diag;   /* Unit or non-unit diagonal */
   ATL_INT M, N;           /* rows (cols) of array B */
   ATL_INT minrhs;         /* # of rhs all jobs get */
   int rb;                 /* block size for RHS division */
   int neblks;             /* # of extra blocks dealt out */
   int rr;                 /* size of last extra-block, or 0 */
};

typedef struct ATL_ttrsm_amm ATL_ttrsm_amm_t;
struct ATL_ttrsm_amm
{
   ammkern_t amm_b0, amm_b1;
   ablk2cmat_t blk2c;
   cm2am_t a2blk, b2blk;
   void *AblkCtr;      /* deals out copy of A blocks column-wise */
   void *rhsCtr;       /* deals out column panels of Right Hand Side (B/X) */
   int *AcpyBV;        /* nablks-len BV; 0: A blk not yet copied */
   void *Acpymut;      /* mutex protecting AcpyBV */
   TYPE *wA;           /* ptr to start of shared A workspace */
   TYPE *w;            /* nthr separate blkszB+panszB workspace */
   const TYPE *A;      /* M-order triangular matrix */
   TYPE *X;            /* MxN matrix of right hand sides (RHS) */
   size_t lda, ldx;
   size_t blkszB;      /* mb*nb, use size_t to prevent overflow with mul */
   size_t blkszA;      /* mb*mb */
   size_t panszC;      /* size of col panel of C : (nmblks-1)*blkszB */
   size_t wsL;         /* # of elts of each core's workspace */
   TYPE alpha;         /* scale to apply to X before access */
   ATL_INT M, N;       /* M: size of A, N: NRHS */
   VINT AcpyDone;      /* set to 1 when copy of entire A complete */
   int nmblks;         /* number of row panels B is split into */
   int nnblks;         /* number of AMM blocks in a row panel */
   int nablks;         /* # of non-diagonal blocks in A */
   int nxblks;         /* # of blocks in X, C has nnblks less */
   int nbf;            /* N%nb, if 0, nb */
   int nnuf;           /* (nbf+nu-1)/nu */
   int mb0;            /* size of 1st diag blk: M%mb, if 0, mb */
   int MB0;            /* for k-vector kernels, (mb0+ku-1)/ku*ku, else mb0 */
   int mb;             /* row & col blocking of A (mb & kb of amm) */
   int nb;             /* column blocking of X/B/C (nb of amm) */
   int mu;             /* M unrolling for amm kernel */
   int nu;             /* N unrolling for amm kernel */
   int ku;             /* K unrolling for amm kernel */
   int nnu, nmu, nmu0;  /* I think can kill nmu0 -> test then KILL */
   enum ATLAS_TRANS TA;
   enum ATLAS_DIAG uplo;
   enum ATLAS_DIAG diag;
};
/*
 * Function prototypes
 */
int Mjoin(PATL,tGetTrsmInfo)(ATL_ttrsm_amm_t *pd, int P, enum ATLAS_TRANS TA,
                             ATL_CSZT M, ATL_CSZT N, const SCALAR beta);
#endif
#endif
