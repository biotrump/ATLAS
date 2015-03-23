#ifndef ATLAS_BCAMM_H
   #define ATLAS_BCAMM_H 1

#include "atlas_misc.h"
#ifdef TYPE
#include "atlas_amm.h"
int Mjoin(PATL,getrf_bcAmm_info)(ATL_INT M, ATL_INT N, int *B, int *R, int *C);
int Mjoin(PATL,tgetrf_bcAMM)
   (ATL_CINT major, ATL_CINT M, ATL_CINT N, TYPE *A, ATL_CINT lda, int *ipiv);
void Mjoin(PATL,bcAblk2cmat)
   (int m, int n, TYPE *W, int nb, TYPE *A, int lda, int nt, ablk2cmat_t ca2cc);
void Mjoin(PATL,bcAm2rm)
   (int m, int n, TYPE *W, int nb, TYPE *A, int lda, int nt, am2cm_t a2rn);
void Mjoin(PATL,bcG2L_cpy)
   (int m, int n, TYPE* A, int lda, TYPE* W, int ldw, int nb, int nt);
void Mjoin(PATL,bcL2G_blkcpy)
   ( int m, int n, TYPE *W, int nb, TYPE *A, int lda, int nt);
void Mjoin(PATL,bcL2G_cpy)
   (int m, int n, TYPE* W, int ldw, TYPE* A, int lda, int nb, int nt);
void Mjoin(PATL,bcRm2am)
   (int m, int n, TYPE *A, int lda, TYPE *W, int nb, int nt, cm2am_t r2an);
#endif
/*
 * This structure is used to encode block-cyclic index info in ipiv using
 * bitfields.  See src/threads/lapack/amm for routines that use it.
 */
typedef struct ATL_bcpiv ATL_bcpiv_t;
struct ATL_bcpiv
{
   ATL_UINT R;       /* # of procs matrix dim distributed over */
   ATL_UINT C;       /* # of procs matrix dim distributed over */
   ATL_UINT B;       /* block factor used to distribute N over P */
   ATL_UINT MU;      /* amm's unrolling along the M dim */
   ATL_UINT NU;      /* amm's unrolling along the N dim */
   ATL_UINT nbRNK;   /* # of bits encoding prank of local block owner */
   ATL_UINT nbSBN;   /* # of bits encoding subblock number */
   ATL_UINT nbSBR;   /* # of bits encoding subblock row index */
   ATL_UINT *ipiv0;  /* ptr to original pivot array */
   ATL_UINT *ipiv;   /* ptr to pivot array being localized */
   ATL_INT inci;     /* piv stride; If inci<0, pivs applied in reverse order */
   ATL_UINT neSB;    /* # of elts in subblock (mu*nu) */
   ATL_UINT neMB;    /* # elts in major block (B*B) */
   void **larrs;     /* ptrs to beginning of local arrays */
   int *lldps;       /* ptrs to leading dimension of local panels */
};
void *ATL_bcIpivInit(int MN, int *ipiv, int inci, int R, int C,
                     int N, int B, int mu, int nu);
void ATL_bcIpivEncode(ATL_bcpiv_t *bp, int n, int I, int iadj);
void ATL_bcIpivDecode(ATL_bcpiv_t *bp, int n, int I);
#define ATL_Free_bcpiv_t(m_) free(m_)

#ifdef TYPE
void Mjoin(PATL,bcLaswp_amm)(ATL_bcpiv_t *bp, ATL_CINT N, ATL_CINT coff,
                             TYPE *A, ATL_CINT lda, ATL_CINT K1, ATL_CINT K2,
                             ATL_CINT lpj, ATL_CINT nnu, ATL_CINT ZOFF);
#endif
#endif
