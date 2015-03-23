/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2014 R. Clint Whaley
 *
 * Code contributers : R. Clint Whaley, Jeff Horner, Antoine Petitet
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
/*
 * =====================================================================
 * Include files
 * =====================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#ifdef GCCWIN
   ___main(){} __main(){} MAIN__(){} _MAIN_(){}
   #ifndef isdigit
      #define isdigit(ch_) ( ((ch_)=='0')||((ch_)=='1')||((ch_)=='2')|| \
                             ((ch_)=='3')||((ch_)=='4')||((ch_)=='5')|| \
                             ((ch_)=='6')||((ch_)=='7')||((ch_)=='8')|| \
                             ((ch_)=='9') )
   #endif
#else
   #include <ctype.h>
#endif

#include "atlas_misc.h"
#include "atlas_tst.h"
/*
 * =====================================================================
 * #define macro constants
 * =====================================================================
 */
#define    MEGA                     1000000.0
#if defined( SREAL ) || defined( SCPLX )
   #define    THRESH                        50.0f
#else
   #define    THRESH                        50.0
#endif
static int NSAMP=1;

/* #define    ATLAS_DEBUG */
/*
 * =====================================================================
 * # macro functions
 * =====================================================================
 *
 * The following and mutually exclusive  macros  allow to select various
 * BLAS implementations to test the ATLAS implementation against:
 *
 *    USE_F77_BLAS     : Fortran 77 BLAS interface,
 *    USE_L3_REFERENCE : C ATLAS reference implementation,
 *
 * If none of these macros is defined at compile time, the  ATLAS imple-
 * mentation is to be tested against itself,  after all this is the only
 * version we are sure to have available.
 *
 * By default the mono-threaded  ATLAS  routines are tested. To test the
 * multi-threaded ATLAS routines, define the following macro:
 *    USE_L3_PTHREADS  : multi-threaded ATLAS implementation.
 */
#define USE_F77_BLAS
#ifdef ATL_USEPTHREADS
#define USE_L3_PTHREADS
#endif

#ifndef USE_L3_PTHREADS
#endif
/*
 * =====================================================================
 */
#if defined(USE_F77_BLAS) || defined(TEST_F77)
   #define  TP3 Mjoin(PATL,f77)
#elif defined( USE_L3_REFERENCE )
   #include "atlas_reflevel3.h"
   #define  TP3      Mjoin( PATL,   ref )
#else /* defined( USE_L3_ATLAS ) */  /* use ATLAS itself !! (default) */
   #include "atlas_level3.h"
   #define  TP3      PATL
#endif

#define trusted_gemm(TA, TB, M, N, K, al, A, lA, B, lB, be, C, lC) \
   Mjoin(TP3, gemm)(TA, TB, M, N, K, al, A, lA, B, lB, be, C, lC)
#define trusted_trsm( SI, UP, TA, DI, M, N, al, A, lA, B, lB) \
   Mjoin( TP3, trsm  )(  SI, UP, TA, DI, M, N,    al, A, lA, B, lB)
#define trusted_trmm(SI, UP, TA, DI, M, N, al, A, lA, B, lB) \
   Mjoin( TP3, trmm  )(SI, UP, TA, DI, M, N,  al, A, lA, B, lB)
#define trusted_syr2k(UP, TA, N, K, al, A, lA, B, lB, be, C, lC) \
   Mjoin( TP3, syr2k )(UP, TA, N, K, al, A, lA, B, lB, be, C, lC)
#define trusted_symm(SI, UP, M, N, al, A, lA, B, lB, be, C, lC) \
   Mjoin(TP3,symm)(SI, UP, M, N, al, A, lA, B, lB, be, C, lC)
#define trusted_syrk(UP, TA, N, K, al, A, lA, be, C, lC) \
   Mjoin(TP3, syrk)(UP, TA, N, K, al, A, lA, be, C, lC)
#ifdef TCPLX
   #define trusted_herk(UP, TA, N, K, al, A, lA, be, C, lC) \
      Mjoin(TP3, herk)(UP, TA, N, K, al, A, lA, be, C, lC)
   #define trusted_her2k(UP, TA, N, K, al, A, lA, B, lB, be, C, lC) \
      Mjoin(TP3,her2k)(UP, TA, N, K, al, A, lA, B, lB, be, C, lC)
   #define trusted_hemm( SI, UP, M, N, al, A, lA, B, lB, be, C, lC) \
      Mjoin(TP3,hemm)(SI, UP, M, N, al, A, lA, B, lB, be, C, lC)
#endif

/*
 * ATLAS version of the BLAS to test.
 */
/* #define PSPRK */
/* #define HPRK */
#if defined(PSPRK) || defined(HPRK)
#include "atlas_pkblas.h"
#endif
#if defined(TEST_F77)
   #define  AP3 Mjoin(PATL,f77)
   #define  AP4 Mjoin(PATL,f77)
#elif defined(USE_L3_PTHREADS)
   #include "atlas_pthreads.h"
   #include "atlas_ptlvl3.h"
   #include "atlas_tlvl3.h"
   #define  AP3  Mjoin(PATL, t)
   #define  AP4      Mjoin( PATL,   t  )
#else
   #include "atlas_level3.h"
   #define  AP3      PATL
   #define  AP4      PATL
#endif

#ifdef THREADS_AvC
   #define test_gemm(TA, TB, M, N, K, al, A, lA, B, lB, be, C, lC) \
      Mjoin(PATL,tgemm)(TA, TB, M, N, K, al, A, lA, B, lB, be, C, lC)
   #define test_trsm(SI, UP, TA, DI, M, N, al, A, lA, B, lB) \
      Mjoin(PATL,ttrsm)(SI, UP, TA, DI, M, N, al, A, lA, B, lB)
   #define test_trmm(SI, UP, TA, DI, M, N, al, A, lA, B, lB) \
      Mjoin(PATL,ttrmm)(SI, UP, TA, DI, M, N, al, A, lA, B, lB)
   #define test_syr2k(UP, TA, N, K, al, A, lA, B, lB, be, C, lC) \
      Mjoin(PATL,tsyr2k)(UP, TA, N, K, al, A, lA, B, lB, be, C, lC)
   #define test_symm(SI, UP, M, N, al, A, lA, B, lB, be, C, lC) \
      Mjoin(PATL,tsymm)(SI, UP, M, N, al, A, lA, B, lB, be, C, lC)
   #ifndef PSPRK
      #define test_syrk(UP, TA, N, K, al, A, lA, be, C, lC) \
         Mjoin(PATL,tsyrk)(UP, TA, N, K, al, A, lA, be, C, lC)
   #else
      #define test_syrk(UP, TA, N, K, al, A, lA, be, C, lC) \
         Mjoin(PATL,sprk)(PackGen, TA, UP, 0, N, K, al, A, 0, 0, lA, \
                          be, C, 0, 0, lC)
   #endif
   #ifdef TCPLX
      #define test_her2k(UP, TA, N, K, al, A, lA, B, lB, be, C, lC) \
         Mjoin(PATL, ther2k)(UP, TA, N, K, al, A, lA, B, lB, be, C, lC)
      #define test_hemm(SI, UP, M, N, al, A, lA, B, lB, be, C, lC) \
         Mjoin(PATL,themm)(SI, UP, M, N, al, A, lA, B, lB, be, C, lC)
      #ifndef HPRK
         #define test_herk(UP, TA, N, K, al, A, lA, be, C, lC) \
            Mjoin(PATL,therk)(UP, TA, N, K, al, A, lA, be, C, lC)
      #else
         #define test_herk(UP, TA, N, K, al, A, lA, be, C, lC) \
            Mjoin(PATL,hprk)(PackGen, TA, UP, 0, N, K, al, A, 0, 0, lA, \
                             be, C, 0, 0, lC)
      #endif
   #endif
#else
   #define test_gemm(TA, TB, M, N, K, al, A, lA, B, lB, be, C, lC) \
      Mjoin( AP4, gemm  )(TA, TB, M, N, K, al, A, lA, B, lB, be, C, lC)
   #define test_trsm(SI, UP, TA, DI, M, N, al, A, lA, B, lB) \
      Mjoin(AP4, trsm)(SI, UP, TA, DI, M, N, al, A, lA, B, lB)
   #define test_trmm(SI, UP, TA, DI, M, N, al, A, lA, B, lB) \
      Mjoin( AP4, trmm  )(SI, UP, TA, DI, M, N, al, A, lA, B, lB)
   #define test_syr2k(UP, TA, N, K, al, A, lA, B, lB, be, C, lC) \
      Mjoin( AP4, syr2k )(UP, TA, N, K, al, A, lA, B, lB, be, C, lC)
   #define test_symm(SI, UP, M, N, al, A, lA, B, lB, be, C, lC) \
      Mjoin(AP4,symm)(SI, UP, M, N, al, A, lA, B, lB, be, C, lC)
   #ifndef PSPRK
      #define test_syrk(UP, TA, N, K, al, A, lA, be, C, lC) \
         Mjoin(AP4,syrk)(UP, TA, N, K, al, A, lA, be, C, lC)
   #else
      #define test_syrk(UP, TA, N, K, al, A, lA, be, C, lC) \
         Mjoin(PATL,sprk)(PackGen, TA, UP, 0, N, K, al, A, 0, 0, lA, \
                          be, C, 0, 0, lC)
   #endif
   #ifdef TCPLX
      #define test_her2k(UP, TA, N, K, al, A, lA, B, lB, be, C, lC) \
         Mjoin( AP4, her2k )(UP, TA, N, K, al, A, lA, B, lB, be, C, lC)
      #define test_hemm(SI, UP, M, N, al, A, lA, B, lB, be, C, lC) \
         Mjoin(AP4,hemm)(SI, UP, M, N, al, A, lA, B, lB, be, C, lC)
      #ifndef HPRK
         #define test_herk(UP, TA, N, K, al, A, lA, be, C, lC) \
            Mjoin(AP4,herk)(UP, TA, N, K, al, A, lA, be, C, lC)
      #else
         #define test_herk(UP, TA, N, K, al, A, lA, be, C, lC) \
            Mjoin(PATL,hprk)(PackGen, TA, UP, 0, N, K, al, A, 0, 0, lA, \
                             be, C, 0, 0, lC)
      #endif
   #endif
#endif

/*
 * =====================================================================
 * macro functions
 * =====================================================================
 */
#ifdef TCPLX
#define Mabs1(X) (Mabs(*X) + Mabs(*(X+1)))
#else
#define Mabs1(X) (Mabs(X))
#endif

#ifdef  ATL_NTHREADS
#define LCSIZE          ATL_NTHREADS * L2SIZE
#else
#define LCSIZE          L2SIZE
#endif
/*
 * =====================================================================
 * typedef definitions
 * =====================================================================
 */
#ifdef TREAL
enum LVL3_ROUT { GEMM=0, SYMM, SYRK, SYR2K, TRMM, TRSM, ALLROUTS };
#else
enum LVL3_ROUT
{ GEMM=0, HEMM, HERK, HER2K, SYMM, SYRK, SYR2K, TRMM, TRSM, ALLROUTS };
#endif
/*
 * =====================================================================
 * Prototypes for the testing routines
 * =====================================================================
 */
double     opbl3
(  const enum LVL3_ROUT,           const int,      const int,
   const int );
void       trddom
(  const enum ATLAS_UPLO,          const int,      TYPE *,
   const int );

TYPE       gemmtst
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const enum ATLAS_TRANS,
   const enum ATLAS_TRANS,         const int,      const int,
   const int,      const SCALAR,   const int,      const int,
   const SCALAR,   const int,      const TYPE,     double *,
   double *,       double *,       double * );
TYPE       symmtst
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const enum ATLAS_SIDE,
   const enum ATLAS_UPLO,          const int,      const int,
   const SCALAR,   const int,      const int,      const SCALAR,
   const int,      const TYPE,     double *,       double *,
   double *,       double * );
TYPE       syr2ktst
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const enum ATLAS_UPLO,
   const enum ATLAS_TRANS,         const int,      const int,
   const SCALAR,   const int,      const int,      const SCALAR,
   const int,      const TYPE,     double *,       double *,
   double *,       double * );
TYPE       syrktst
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const enum ATLAS_UPLO,
   const enum ATLAS_TRANS,         const int,      const int,
   const SCALAR,   const int,      const SCALAR,   const int,
   const TYPE,     double *,       double *,       double *,
   double * );
TYPE       trmmtst
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const enum ATLAS_SIDE,
   const enum ATLAS_UPLO,          const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,          const int,      const int,
   const SCALAR,   const int,      const int,      const TYPE,
   double *,       double *,       double *,       double * );
TYPE       trsmtst
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const enum ATLAS_SIDE,
   const enum ATLAS_UPLO,          const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,          const int,      const int,
   const SCALAR,   const int,      const int,      const TYPE,
   double *,       double *,       double *,       double * );

int        gemmcase
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const int,
   const enum ATLAS_TRANS,         const enum ATLAS_TRANS,
   const int,      const int,      const int,      const SCALAR,
   const int,      const int,      const SCALAR,   const int,
   const TYPE,     double *,       double *,       double *,
   double * );
int        symmcase
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const int,
   const enum ATLAS_SIDE,          const enum ATLAS_UPLO,
   const int,      const int,      const SCALAR,   const int,
   const int,      const SCALAR,   const int,      const TYPE,
   double *,       double *,       double *,       double * );
int        syr2kcase
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const int,
   const enum ATLAS_UPLO,          const enum ATLAS_TRANS,
   const int,      const int,      const SCALAR,   const int,
   const int,      const SCALAR,   const int,      const TYPE,
   double *,       double *,       double *,       double * );
int        syrkcase
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const int,
   const enum ATLAS_UPLO,          const enum ATLAS_TRANS,
   const int,      const int,      const SCALAR,   const int,
   const SCALAR,   const int,      const TYPE,     double *,
   double *,       double *,       double * );
int        trxmcase
(  const int CACHESIZE,
   const enum LVL3_ROUT,           const int,      const int,
   const enum ATLAS_SIDE,          const enum ATLAS_UPLO,
   const enum ATLAS_TRANS,         const enum ATLAS_DIAG,
   const int,      const int,      const SCALAR,   const int,
   const int,      const TYPE,     double *,       double *,
   double *,       double * );


int        main
(  int,            char ** );
/*
 * =====================================================================
 */
double opbl3
(
   const enum LVL3_ROUT       ROUT,
   const int                  M,
   const int                  N,
   const int                  K
)
{
   double                     adds = 0.0, em, en, ek, muls = 0.0;
/*
 * On entry,  M,  N,  and K contain parameter values used by the Level 3
 * BLAS.  The output matrix is always M x N or N x N if symmetric, but K
 * has different uses in different contexts. For example, in the matrix-
 * matrix multiply routine,  we  have C = A * B where  C is M x N,  A is
 * M x K, and B is K x N. In xSYMM, xHEMM, xTRMM, and xTRSM, K indicates
 * whether the matrix A is applied on the left or right.  If K <= 0, the
 * matrix is aqpplied on the left, and if K > 0, on  the  right.
 */
   if( M <= 0 ) return( 0.0 );

   em = (double)(M); en = (double)(N); ek = (double)(K);

   if(      ROUT == GEMM ) { muls = em * ek * en; adds = em * ek * en; }
#ifdef TREAL
   else if( ROUT == SYMM )
#else
   else if( ( ROUT == SYMM ) || ( ROUT == HEMM ) )
#endif
   {                 /* If K <= 0, assume A multiplies B on the left. */
      if( K <= 0 ) { muls = em * em * en; adds = em * em * en; }
      else         { muls = em * en * en; adds = em * en * en; }
   }
   else if( ROUT == TRMM )
   {                 /* If K <= 0, assume A multiplies B on the left. */
      if( K <= 0 )
      {
         muls = en * em * ( em + 1.0 ) / 2.0;
         adds = en * em * ( em - 1.0 ) / 2.0;
      }
      else
      {
         muls = em * en * ( en + 1.0 ) / 2.0;
         adds = em * en * ( en - 1.0 ) / 2.0;
      }
   }
#ifdef TREAL
   else if( ROUT == SYRK )
#else
   else if( ( ROUT == SYRK ) || ( ROUT == HERK ) )
#endif
   {
      muls = ek * em * ( em + 1.0 ) / 2.0;
      adds = ek * em * ( em + 1.0 ) / 2.0;
   }
#ifdef TREAL
   else if( ROUT == SYR2K )
#else
   else if( ( ROUT == SYR2K ) || ( ROUT == HER2K ) )
#endif
   { muls = ek * em * em; adds = ek * em * em + em; }
   else if( ROUT == TRSM )
   {                 /* If K <= 0, assume A multiplies B on the left. */
      if( K <= 0 )
      {
         muls = en * em * ( em + 1.0 ) / 2.0;
         adds = en * em * ( em - 1.0 ) / 2.0;
      }
      else
      {
         muls = em * en * ( en + 1.0 ) / 2.0;
         adds = em * en * ( en - 1.0 ) / 2.0;
      }
   }
#ifdef TREAL
   return(       muls +       adds );
#else
   return( 6.0 * muls + 2.0 * adds );
#endif
}

void trddom
(
   const enum ATLAS_UPLO      UPLO,
   const int                  N,
   TYPE                       * A,
   const int                  LDA
)
{
/*
 * Scale strictly lower (resp. upper) part of triangular matrix by 1 / N
 * to make it diagonally dominant.
 */
   int                        i, iaij, j, jaj, lda2 = ( LDA SHIFT ),
                              ldap12 = (( LDA + 1 ) SHIFT);
   TYPE                       alpha;

   if( N <= 0 ) return;

   alpha = ATL_rone / (TYPE)(N);

   if( UPLO == AtlasUpper )
   {
      for( j = 0, jaj = 0; j < N; j++, jaj += lda2 )
      {
         for( i = 0, iaij = jaj; i < j; i++, iaij += (1 SHIFT) )
         {
            A[iaij  ] *= alpha;
#ifdef TCPLX
            A[iaij+1] *= alpha;
#endif
         }
         if( A[iaij  ] >= ATL_rzero ) A[iaij  ] += ATL_rone;
         else                         A[iaij  ] -= ATL_rone;
#ifdef TCPLX
         if( A[iaij+1] >= ATL_rzero ) A[iaij+1] += ATL_rone;
         else                         A[iaij+1] -= ATL_rone;
#endif
      }
   }
   else
   {
      for( j = N-1, jaj = (N-1)*ldap12; j >= 0; j--, jaj -= ldap12 )
      {
         if( A[jaj  ] >= ATL_rzero ) A[jaj  ] += ATL_rone;
         else                        A[jaj  ] -= ATL_rone;
#ifdef TCPLX
         if( A[jaj+1] >= ATL_rzero ) A[jaj+1] += ATL_rone;
         else                        A[jaj+1] -= ATL_rone;
#endif
         for( i = j+1, iaij = jaj+(1 SHIFT); i < N; i++, iaij += (1 SHIFT) )
         {
            A[iaij  ] *= alpha;
#ifdef TCPLX
            A[iaij+1] *= alpha;
#endif
         }
      }
   }
}
/*
 * =====================================================================
 * tst functions
 * =====================================================================
 */
TYPE gemmtst
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const enum ATLAS_TRANS     TRANSA,
   const enum ATLAS_TRANS     TRANSB,
   const int                  M,
   const int                  N,
   const int                  K,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const SCALAR               BETA,
   const int                  LDC,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     l2ret, ops, t0, ttest, ttrust;
   TYPE                       normA, normB, normC, normD, resid;
   TYPE                       * A  = NULL, * B = NULL, * C = NULL, * C0,
                              * a, * b, * c;
   int                        mA, mB, nA, nB, Aseed, Bseed, Cseed;
   const int DOFLUSH=(ATL_MulBySize((size_t)M*N + (size_t)K*(M+N)) <
                      ((size_t)CACHESIZE)*4);

   *TTRUST0 = *TTEST0 = *MFTEST0 = *MFTRUST0 = 0.0;
   if( ( M == 0 ) || ( N == 0 ) ) { return( ATL_rzero ); }

   if( TRANSA == AtlasNoTrans ) { mA = M; nA = K; }
   else                         { mA = K; nA = M; }
   if( TRANSB == AtlasNoTrans ) { mB = K; nB = N; }
   else                         { mB = N; nB = K; }

   ops = opbl3( ROUT, M, N, K );
/*
 * Allocate L2 cache space, A, X, Y and Y0
 */
   if (DOFLUSH)
      l2ret = ATL_flushcache( CACHESIZE );
   A  = (TYPE *)malloc( ATL_MulBySize( LDA ) * nA     );
   B  = (TYPE *)malloc( ATL_MulBySize( LDB ) * nB     );
   C  = (TYPE *)malloc( ATL_MulBySize( LDC ) * N  * 2 );

   if( ( A == NULL ) || ( B == NULL ) || ( C == NULL ) )
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( 0 );
      if( A ) free( A );
      if( B ) free( B );
      if( C ) free( C );
      return( ATL_rnone );
   }

   C0 = C + LDC * ( N SHIFT );
/*
 * Generate random operands
 */
   Aseed = mA * nA + 513 *  7 + 90;
   Bseed = mB * nB + 127 * 50 + 77;
   Cseed = M  * N  + 101 *  2 + 53;

   Mjoin( PATL, gegen )( mA, nA, A,  LDA, Aseed );
   Mjoin( PATL, gegen )( mB, nB, B,  LDB, Bseed );
   Mjoin( PATL, gegen )( M,  N,  C,  LDC, Cseed );
   Mjoin( PATL, gegen )( M,  N,  C0, LDC, Cseed );
/*
 * Compute the norm of C for later use in testing
 */
   if( TEST )
   {
      normC = Mjoin( PATL, genrm1 )( M, N, C, LDC );
      if( Mabs1( BETA ) > ATL_rone ) normC *= Mabs1( BETA  );
      if( normC == ATL_rzero ) normC = ATL_rone;
   }
   else { normC = ATL_rone; }
/*
 * Start cold cache timing operations for the trusted routine
 */
   a = A; b = B; c = C0;

   if (DOFLUSH)
      l2ret  = ATL_flushcache( -1 );
   t0     = time00();
   trusted_gemm( TRANSA, TRANSB, M, N, K, ALPHA, a, LDA, b, LDB, BETA, c, LDC );
   ttrust = time00() - t0;
   if( ttrust > 0.0 )
   { *TTRUST0 = ttrust; *MFTRUST0 = ops / ( ttrust * MEGA ); }
   if (DOFLUSH)
      l2ret  = ATL_flushcache( 0 );
   free(A);
   free(B);
   free(C);
   return( ATL_rzero );
}

TYPE symmtst
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const SCALAR               BETA,
   const int                  LDC,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     l2ret, ops, t0, ttest, ttrust;
   TYPE                       normA, normB, normC, normD, resid;
   TYPE                       * A  = NULL, * B = NULL, * C = NULL, * C0,
                              * a, * b, * c;
   int                        nA, Aseed, Bseed, Cseed;
   const size_t szA = (SIDE == AtlasLeft) ? M : N;
   const int DOFLUSH=(ATL_MulBySize((size_t)M*N*2 + szA)<((size_t)CACHESIZE)*4);

   *TTRUST0 = *TTEST0 = *MFTEST0 = *MFTRUST0 = 0.0;
   if( N == 0 ) { return( ATL_rzero ); }

   if( SIDE == AtlasLeft ) { ops = opbl3( ROUT, M, N, -1 ); nA = M; }
   else                    { ops = opbl3( ROUT, M, N,  1 ); nA = N; }
/*
 * Allocate L2 cache space, A, X, Y and Y0
 */
   if (DOFLUSH)
      l2ret = ATL_flushcache( CACHESIZE );
   A  = (TYPE *)malloc( ATL_MulBySize( LDA ) * nA     );
   B  = (TYPE *)malloc( ATL_MulBySize( LDB ) * N      );
   C  = (TYPE *)malloc( ATL_MulBySize( LDC ) * N  * 2 );

   if( ( A == NULL ) || ( B == NULL ) || ( C == NULL ) )
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( 0 );
      if( A ) free( A );
      if( B ) free( B );
      if( C ) free( C );
      return( ATL_rnone );
   }

   C0 = C + LDC * ( N SHIFT );
/*
 * Generate random operands
 */
   Aseed = nA * nA + 513 *  7 + 90;
   Bseed = M  * N  + 127 * 50 + 77;
   Cseed = M  * N  + 101 *  2 + 53;

   Mjoin( PATL, gegen )( nA, nA, A,  LDA, Aseed );
   Mjoin( PATL, gegen )( M,  N,  B,  LDB, Bseed );
   Mjoin( PATL, gegen )( M,  N,  C,  LDC, Cseed );
   Mjoin( PATL, gegen )( M,  N,  C0, LDC, Cseed );
/*
 * Compute the norm of C for later use in testing
 */
   if( TEST )
   {
      normC = Mjoin( PATL, genrm1 )( M, N, C, LDC );
      if( Mabs1( BETA ) > ATL_rone ) normC *= Mabs1( BETA  );
      if( normC == ATL_rzero ) normC = ATL_rone;
   }
   else { normC = ATL_rone; }
/*
 * Start cold cache timing operations for the trusted routine
 */
   a = A; b = B; c = C0;

#ifdef TREAL
   if (DOFLUSH)
      l2ret = ATL_flushcache( -1 );
   t0     = time00();
   trusted_symm( SIDE, UPLO, M, N, ALPHA, a, LDA, b, LDB, BETA, c, LDC );
   ttrust = time00() - t0;
   if( ttrust > 0.0 )
   { *TTRUST0 = ttrust; *MFTRUST0 = ops / ( ttrust * MEGA ); }
#else
   if( ROUT == SYMM )
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( -1 );
      t0     = time00();
      trusted_symm( SIDE, UPLO, M, N, ALPHA, a, LDA, b, LDB, BETA, c, LDC );
      ttrust = time00() - t0;
      if( ttrust > 0.0 )
      { *TTRUST0 = ttrust; *MFTRUST0 = ops / ( ttrust * MEGA ); }
   }
   else /* if( ROUT == HEMM ) */
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( -1 );
      t0     = time00();
      trusted_hemm( SIDE, UPLO, M, N, ALPHA, a, LDA, b, LDB, BETA, c, LDC );
      ttrust = time00() - t0;
      if( ttrust > 0.0 )
      { *TTRUST0 = ttrust; *MFTRUST0 = ops / ( ttrust * MEGA ); }
   }
#endif
   if (DOFLUSH)
      l2ret  = ATL_flushcache( 0 );
   free(A);
   free(B);
   free(C);
   return( ATL_rzero );
}

TYPE syr2ktst
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const int                  N,
   const int                  K,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const SCALAR               BETA,
   const int                  LDC,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     l2ret, ops, t0, ttest, ttrust;
   TYPE                       normA, normB, normC, normD, resid;
   TYPE                       * A = NULL, * B = NULL, * C = NULL, * C0,
                              * a, * b, * c;
   int                        mAB, nAB, Aseed, Bseed, Cseed;
   enum ATLAS_TRANS           ta;
   const int DOFLUSH=(ATL_MulBySize((size_t)N*N + (size_t)N*K*2)
                      < ((size_t)CACHESIZE)*4);

   *TTRUST0 = *TTEST0 = *MFTEST0 = *MFTRUST0 = 0.0;
   if( N == 0 ) { return( ATL_rzero ); }

   if( TRANS == AtlasNoTrans )
   { ta = TRANS; mAB = N; nAB = K; }
   else
   { ta = ( ROUT == SYR2K ? AtlasTrans : AtlasConjTrans ); mAB = K; nAB = N; }

   ops = opbl3( ROUT, N, 0, K );
/*
 * Allocate L2 cache space, A, C and C0
 */
   if (DOFLUSH)
      l2ret = ATL_flushcache( CACHESIZE );
   A = (TYPE *)malloc( ATL_MulBySize( LDA ) * nAB    );
   B = (TYPE *)malloc( ATL_MulBySize( LDB ) * nAB    );
   C = (TYPE *)malloc( ATL_MulBySize( LDC ) * N  * 2 );

   if( ( A == NULL ) || ( B == NULL ) || ( C == NULL ) )
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( 0 );
      if( A ) free( A );
      if( B ) free( B );
      if( C ) free( C );
      return( ATL_rnone );
   }

   C0 = C + LDC * ( N SHIFT );
/*
 * Generate random operands
 */
   Aseed = mAB * nAB + 513 *  7 + 90;
   Bseed = mAB * nAB + 127 * 50 + 77;
   Cseed = N   * N   + 101 *  2 + 53;

   Mjoin( PATL, gegen )( mAB, nAB, A,  LDA, Aseed );
   Mjoin( PATL, gegen )( mAB, nAB, B,  LDB, Bseed );
   Mjoin( PATL, gegen )( N,   N,   C,  LDC, Cseed );
   Mjoin( PATL, gegen )( N,   N,   C0, LDC, Cseed );
/*
 * Compute the norm of C for later use in testing
 */
   if( TEST )
   {
#ifdef TREAL
      normC = Mjoin( PATL, synrm )( UPLO, N, C, LDC );
      if( Mabs1( BETA ) > ATL_rone ) normC *= Mabs1( BETA  );
      if( normC == ATL_rzero ) normC = ATL_rone;
#else
      if( ROUT == SYR2K )
      {
         normC = Mjoin( PATL, synrm )( UPLO, N, C, LDC );
         if( Mabs1( BETA ) > ATL_rone ) normC *= Mabs1( BETA  );
         if( normC == ATL_rzero ) normC = ATL_rone;
      }
      else
      {
         normC = Mjoin( PATL, henrm )( UPLO, N, C, LDC );
         if( Mabs( BETA[0] ) > ATL_rone ) normC *= Mabs( BETA[0] );
         if( normC == ATL_rzero ) normC = ATL_rone;
      }
#endif
   }
   else { normC = ATL_rone; }
/*
 * Start cold cache timing operations for the trusted routine
 */
   a = A; b = B; c = C0;
#ifdef TREAL
   if (DOFLUSH)
      l2ret = ATL_flushcache( -1 );
   t0     = time00();
   trusted_syr2k( UPLO, ta, N, K, ALPHA, a, LDA, b, LDB, BETA, c, LDC );
   ttrust = time00() - t0;
#else
   if( ROUT == SYR2K )
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( -1 );
      t0     = time00();
      trusted_syr2k( UPLO, ta, N, K, ALPHA, a, LDA, b, LDB, BETA, c, LDC );
      ttrust = time00() - t0;
   }
   else /* if( ROUT == HER2K ) */
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( -1 );
      t0     = time00();
      trusted_her2k( UPLO, ta, N, K, ALPHA, a, LDA, b, LDB, (TYPE)(BETA[0]),
                     c, LDC );
      ttrust = time00() - t0;
   }
#endif
   if( ttrust > 0.0 )
   { *TTRUST0 = ttrust; *MFTRUST0 = ops / ( ttrust * MEGA ); }
   if (DOFLUSH)
      l2ret  = ATL_flushcache( 0 );
   free(A);
   free(B);
   free(C);
   return( ATL_rzero );
}

TYPE syrktst
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const int                  N,
   const int                  K,
   const SCALAR               ALPHA,
   const int                  LDA,
   const SCALAR               BETA,
   const int                  LDC,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     l2ret, ops, t0, ttest, ttrust;
   TYPE                       normA, normC, normD, resid;
   TYPE                       * A = NULL, * C = NULL, * C0, * a, * c;
   int                        mA, nA, Aseed, Cseed;
   enum ATLAS_TRANS           ta;
   const int DOFLUSH=(ATL_MulBySize((size_t)N*N + (size_t)N*K)
                      < ((size_t)CACHESIZE)*4);

   *TTRUST0 = *TTEST0 = *MFTEST0 = *MFTRUST0 = 0.0;
   if( N == 0 ) { return( ATL_rzero ); }

   if( TRANS == AtlasNoTrans )
   { ta = TRANS; mA = N; nA = K; }
   else
   { ta = ( ROUT == SYRK ? AtlasTrans : AtlasConjTrans ); mA = K; nA = N; }

   ops = opbl3( ROUT, N, 0, K );
/*
 * Allocate L2 cache space, A, C and C0
 */
   if (DOFLUSH)
      l2ret = ATL_flushcache( CACHESIZE );
   A = (TYPE *)malloc( ATL_MulBySize( LDA ) * nA     );
   C = (TYPE *)malloc( ATL_MulBySize( LDC ) * N  * 2 );

   if( ( A == NULL ) || ( C == NULL ) )
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( 0 );
      if( A ) free( A );
      if( C ) free( C );
      return( ATL_rnone );
   }

   C0 = C + LDC * ( N SHIFT );
/*
 * Generate random operands
 */
   Aseed = mA * nA + 513 *  7 + 90;
   Cseed = N  * N  + 101 *  2 + 53;

   Mjoin( PATL, gegen )( mA, nA, A,  LDA, Aseed );
   Mjoin( PATL, gegen )( N,  N,  C,  LDC, Cseed );
   Mjoin( PATL, gegen )( N,  N,  C0, LDC, Cseed );
/*
 * Compute the norm of C for later use in testing
 */
   if( TEST )
   {
#ifdef TREAL
      normC = Mjoin( PATL, synrm )( UPLO, N, C, LDC );
      if( Mabs1( BETA ) > ATL_rone ) normC *= Mabs1( BETA  );
      if( normC == ATL_rzero ) normC = ATL_rone;
#else
      if( ROUT == SYRK )
      {
         normC = Mjoin( PATL, synrm )( UPLO, N, C, LDC );
         if( Mabs1( BETA ) > ATL_rone ) normC *= Mabs1( BETA  );
         if( normC == ATL_rzero ) normC = ATL_rone;
      }
      else
      {
         normC = Mjoin( PATL, henrm )( UPLO, N, C, LDC );
         if( Mabs( BETA[0] ) > ATL_rone ) normC *= Mabs( BETA[0] );
         if( normC == ATL_rzero ) normC = ATL_rone;
      }
#endif
   }
   else { normC = ATL_rone; }
/*
 * Start cold cache timing operations for the trusted routine
 */
   a = A; c = C0;
#ifdef TREAL
   if (DOFLUSH)
      l2ret  = ATL_flushcache( -1 );
   t0     = time00();
   trusted_syrk( UPLO, ta, N, K, ALPHA, a, LDA, BETA, c, LDC );
   ttrust = time00() - t0;
#else
   if( ROUT == SYRK )
   {
      if (DOFLUSH)
         l2ret  = ATL_flushcache( -1 );
      t0     = time00();
      trusted_syrk( UPLO, ta, N, K, ALPHA, a, LDA, BETA, c, LDC );
      ttrust = time00() - t0;
   }
   else /* if( ROUT == HERK ) */
   {
      if (DOFLUSH)
         l2ret  = ATL_flushcache( -1 );
      t0     = time00();
      trusted_herk( UPLO, ta, N, K, (TYPE)(ALPHA[0]), a, LDA, (TYPE)(BETA[0]),
                    c, LDC );
      ttrust = time00() - t0;
   }
#endif
   if( ttrust > 0.0 )
   { *TTRUST0 = ttrust; *MFTRUST0 = ops / ( ttrust * MEGA ); }
   if (DOFLUSH)
      l2ret  = ATL_flushcache( 0 );
   free(A);
   free(C);
   return( ATL_rzero );
}

TYPE trmmtst
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_DIAG      DIAG,
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     l2ret, ops, t0, ttest, ttrust;
   TYPE                       normA, normB, normD, resid;
   TYPE                       * A = NULL, * B = NULL, * B0, * a, * b;
   int                        nA, Aseed, Bseed;
   const size_t szA = (SIDE == AtlasLeft) ? M : N;
   const int DOFLUSH=(ATL_MulBySize((size_t)M*N + (size_t)szA*szA)
                      < ((size_t)CACHESIZE)*4);

   *TTRUST0 = *TTEST0 = *MFTEST0 = *MFTRUST0 = 0.0;
   if( ( M == 0 ) || ( N == 0 ) ) { return( ATL_rzero ); }

   if( SIDE == AtlasLeft ) { nA = M; ops = opbl3( ROUT, M, N, -1 ); }
   else                    { nA = N; ops = opbl3( ROUT, M, N,  1 ); }
/*
 * Allocate L2 cache space, A, X and X0
 */
   if (DOFLUSH)
      l2ret = ATL_flushcache( CACHESIZE );
   A = (TYPE *)malloc( ATL_MulBySize( LDA ) * nA    );
   B = (TYPE *)malloc( ATL_MulBySize( LDB ) * N * 2 );

   if( ( A == NULL ) || ( B == NULL ) )
   {
      if (DOFLUSH)
         l2ret  = ATL_flushcache( 0 );
      if( A ) free( A );
      if( B ) free( B );
      return( ATL_rnone );
   }

   B0 = B + LDB * ( N SHIFT );
/*
 * Generate random operands
 */
   Aseed = nA * nA + 513 *  7 + 90;
   Bseed = M  * N  + 127 * 50 + 77;

   Mjoin( PATL, gegen )( nA, nA, A,  LDA, Aseed );
   Mjoin( PATL, gegen )( M,  N,  B,  LDB, Bseed );
   Mjoin( PATL, gegen )( M,  N,  B0, LDB, Bseed );
/*
 * Compute the norm of B for later use in testing
 */
   if( TEST )
   {
      normB = Mjoin( PATL, genrm1 )( M, N, B, LDB );
      if( Mabs1( ALPHA ) > ATL_rone ) normB *= Mabs1( ALPHA );
      if( normB == ATL_rzero ) normB = ATL_rone;
   }
   else { normB = ATL_rone; }
/*
 * Start cold cache timing operations for the trusted routine
 */
   a = A; b = B0;

   if (DOFLUSH)
      l2ret  = ATL_flushcache( -1 );
   t0     = time00();
   trusted_trmm( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, a, LDA, b, LDB );
   ttrust = time00() - t0;
   if( ttrust > 0.0 )
   { *TTRUST0 = ttrust; *MFTRUST0 = ops / ( ttrust * MEGA ); }
   if (DOFLUSH)
      l2ret  = ATL_flushcache( 0 );
   free(A);
   free(B);
   return( ATL_rzero );
}

TYPE trsmtst
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_DIAG      DIAG,
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     l2ret, ops, t0, ttest, ttrust;
   TYPE                       normA, normB, normD, resid;
   TYPE                       * A = NULL, * B = NULL, * B0, * a, * b;
   int                        nA, Aseed, Bseed;
   const size_t szA = (SIDE == AtlasLeft) ? M : N;
   const int DOFLUSH=(ATL_MulBySize((size_t)M*N + (size_t)szA*szA)
                      < ((size_t)CACHESIZE)*4);

   *TTRUST0 = *TTEST0 = *MFTEST0 = *MFTRUST0 = 0.0;
   if( ( M == 0 ) || ( N == 0 ) ) { return( ATL_rzero ); }

   if( SIDE == AtlasLeft ) { nA = M; ops = opbl3( ROUT, M, N, -1 ); }
   else                    { nA = N; ops = opbl3( ROUT, M, N,  1 ); }
/*
 * Allocate L2 cache space, A, X and X0
 */
   if (DOFLUSH)
      l2ret = ATL_flushcache( CACHESIZE );
   A  = (TYPE *)malloc( ATL_MulBySize( LDA ) * nA    );
   B  = (TYPE *)malloc( ATL_MulBySize( LDB ) * N * 2 );

   if( ( A == NULL ) || ( B == NULL ) )
   {
      if (DOFLUSH)
         l2ret = ATL_flushcache( 0 );
      if( A ) free( A );
      if( B ) free( B );
      return( ATL_rnone );
   }

   B0 = B + LDB * ( N SHIFT );
/*
 * Generate random operands
 */
   Aseed = nA * nA + 513 *  7 + 90;
   Bseed = M  * N  + 127 * 50 + 77;

   Mjoin( PATL, gegen )( nA, nA, A,  LDA, Aseed ); trddom( UPLO, nA, A, LDA );
   Mjoin( PATL, gegen )( M,  N,  B,  LDB, Bseed );
   Mjoin( PATL, gegen )( M,  N,  B0, LDB, Bseed );
/*
 * Compute the norm of B for later use in testing
 */
   if( TEST )
   {
      normB = Mjoin( PATL, genrm1 )( M, N, B, LDB );
      if( Mabs1( ALPHA ) > ATL_rone ) normB *= Mabs1( ALPHA );
      if( normB == ATL_rzero ) normB = ATL_rone;
   }
   else { normB = ATL_rone; }
/*
 * Start cold cache timing operations for the trusted routine
 */
   a = A; b = B0;

   if (DOFLUSH)
      l2ret  = ATL_flushcache( -1 );
   t0     = time00();
   trusted_trsm( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, a, LDA, b, LDB );
   ttrust = time00() - t0;
   if( ttrust > 0.0 )
   { *TTRUST0 = ttrust; *MFTRUST0 = ops / ( ttrust * MEGA ); }
   if (DOFLUSH)
      l2ret  = ATL_flushcache( 0 );
   free(A);
   free(B);
   return( ATL_rzero );
}
/*
 * =====================================================================
 * case functions
 * =====================================================================
 */
int gemmcase
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const int                  MFLOP,
   const enum ATLAS_TRANS     TRANSA,
   const enum ATLAS_TRANS     TRANSB,
   const int                  M,
   const int                  N,
   const int                  K,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const SCALAR               BETA,
   const int                  LDC,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     flops, ttrust, ttest, mftrust, mftest, t0;
   TYPE                       resid = ATL_rzero;
#ifdef TREAL
   TYPE                       bet,  beta,    nbeta;
#else
   TYPE                       *bet, beta[2], nbeta[2];
#endif
   TYPE                       * a, * stA, *b, * stB, * c, * stC, * A,
                              * A0 = NULL, * B, * B0 = NULL, * C, * C0 = NULL;
   unsigned long              ir, reps;
   int                        inca, incb, incc, lA, lB, lC, mA, nA, mB, nB,
                              passed, Aseed, Bseed, Cseed;

   if( ( MEGA * MFLOP <= ( flops = opbl3( ROUT, M, N, K ) ) ) || ( TEST ) )
   {
      resid = gemmtst( CACHESIZE, ROUT, TEST, TRANSA, TRANSB, M, N, K, ALPHA,
		       LDA, LDB, BETA, LDC, EPSILON, TTRUST0, TTEST0,
		       MFTRUST0, MFTEST0 );
      if( resid > THRESH ) (void) fprintf( stderr, "   resid=%f\n", resid );
   }
   if( resid < ATL_rzero ) passed = -1;
   else                    passed = ( resid < THRESH );

   if( MEGA * MFLOP <= flops ) return( passed );

   if( TRANSA == AtlasNoTrans ) { mA = M; nA = K; } else { mA = K; nA = M; }
   if( TRANSB == AtlasNoTrans ) { mB = K; nB = N; } else { mB = N; nB = K; }

   inca = LDA  * ( nA SHIFT );
   incb = LDB  * ( nB SHIFT );
   incc = LDC  * ( N  SHIFT );

   lA = inca  * ( ( ATL_DivBySize( LCSIZE ) + mA*nA - 1 ) / ( mA * nA ) );
   lB = incb  * ( ( ATL_DivBySize( LCSIZE ) + mB*nB - 1 ) / ( mB * nB ) );
   lC = incc  * ( ( ATL_DivBySize( LCSIZE ) + M * N - 1 ) / ( M  * N  ) );

   A0 = (TYPE *)malloc( ATL_MulBySize( lA ) );
   B0 = (TYPE *)malloc( ATL_MulBySize( lB ) );
   C0 = (TYPE *)malloc( ATL_MulBySize( lC ) );

   if( ( A0 == NULL ) || ( B0 == NULL ) || ( C0 == NULL ) )
   {
      if( A0 ) free( A0 );
      if( B0 ) free( B0 );
      if( C0 ) free( C0 );
      return( -1 );
   }

   A = A0; stA = A0 + ( lA SHIFT );
   B = B0; stB = B0 + ( lB SHIFT );
   C = C0; stC = C0 + ( lC SHIFT );

#ifdef TREAL
   beta   =  BETA;
   nbeta  = -BETA;
#else
   *beta  =    *BETA; beta [1] =  BETA[1];
   *nbeta = -(*BETA); nbeta[1] = -BETA[1];
#endif

   Aseed = mA * nA + 513 *  7 + 90;
   Bseed = mB * nB + 127 * 50 + 77;
   Cseed = M  * N  + 101 *  2 + 53;

   reps  = ( MEGA * MFLOP ) / flops;
/*
 * Generate the random data and time the trusted routine
 */
   bet = beta; a = A; b = B; c = C;

   Mjoin( PATL, gegen )( lA, 1, A0, lA, Aseed );
   Mjoin( PATL, gegen )( lB, 1, B0, lB, Bseed );
   Mjoin( PATL, gegen )( lC, 1, C0, lC, Cseed );

   t0 = time00();
   for( ir = reps; ir; ir-- )
   {
      trusted_gemm( TRANSA, TRANSB, M, N, K, ALPHA, a, LDA, b, LDB,
                    (SCALAR)(bet), c, LDC );
      a += inca; if( a == stA ) { a = A; }
      b += incb; if( b == stB ) { b = B; }
      c += incc;
      if( c == stC ) { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
   }
   ttrust = time00() - t0;
   if( ttrust > 0.0 ) mftrust = ( reps * flops ) / ( MEGA * ttrust );
   else               mftrust = 0.0;
   ttrust /= reps; *TTRUST0 = ttrust; *MFTRUST0 = mftrust;
   free(A0);
   free(B0);
   free(C0);
   return(0);
}

int symmcase
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const int                  MFLOP,
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const SCALAR               BETA,
   const int                  LDC,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     flops, ttrust, ttest, mftrust, mftest, t0;
   TYPE                       resid = ATL_rzero;
#ifdef TREAL
   TYPE                       bet,  beta,    nbeta;
#else
   TYPE                       *bet, beta[2], nbeta[2];
#endif
   TYPE                       * a, * stA, *b, * stB, * c, * stC, * A,
                              * A0 = NULL, * B, * B0 = NULL, * C, * C0 = NULL;
   unsigned long              ir, reps;
   int                        inca, incb, incc, lA, lB, lC, nA, passed, Aseed,
                              Bseed, Cseed;

   flops = opbl3( ROUT, M, N, ( SIDE == AtlasLeft ? -1 : 1 ) );

   if( ( MEGA * MFLOP <= flops ) || ( TEST ) )
   {
      resid = symmtst( CACHESIZE, ROUT, TEST, SIDE, UPLO, M, N, ALPHA,
		       LDA, LDB, BETA, LDC, EPSILON, TTRUST0, TTEST0,
		       MFTRUST0, MFTEST0 );
      if( resid > THRESH ) (void) fprintf( stderr, "   resid=%f\n", resid );
   }
   if( resid < ATL_rzero ) passed = -1;
   else                    passed = ( resid < THRESH );

   if( MEGA * MFLOP <= flops ) return( passed );

   if( SIDE == AtlasLeft ) { nA = M; } else { nA = N; }

   inca = LDA  * ( nA SHIFT );
   incb = LDB  * ( N  SHIFT );
   incc = LDC  * ( N  SHIFT );

   lA = inca  * ( ( ATL_DivBySize( LCSIZE ) + nA*nA - 1 ) / ( nA * nA ) );
   lB = incb  * ( ( ATL_DivBySize( LCSIZE ) + M * N - 1 ) / ( M  * N  ) );
   lC = incc  * ( ( ATL_DivBySize( LCSIZE ) + M * N - 1 ) / ( M  * N  ) );

   A0 = (TYPE *)malloc( ATL_MulBySize( lA ) );
   B0 = (TYPE *)malloc( ATL_MulBySize( lB ) );
   C0 = (TYPE *)malloc( ATL_MulBySize( lC ) );

   if( ( A0 == NULL ) || ( B0 == NULL ) || ( C0 == NULL ) )
   {
      if( A0 ) free( A0 );
      if( B0 ) free( B0 );
      if( C0 ) free( C0 );
      return( -1 );
   }

   A = A0; stA = A0 + ( lA SHIFT );
   B = B0; stB = B0 + ( lB SHIFT );
   C = C0; stC = C0 + ( lC SHIFT );

#ifdef TREAL
   beta   =  BETA;
   nbeta  = -BETA;
#else
   *beta  =    *BETA; beta [1] =  BETA[1];
   *nbeta = -(*BETA); nbeta[1] = -BETA[1];
#endif

   Aseed = nA * nA + 513 *  7 + 90;
   Bseed = M  * N  + 127 * 50 + 77;
   Cseed = M  * N  + 101 *  2 + 53;

   reps  = ( MEGA * MFLOP ) / flops;
/*
 * Generate the random data and time the trusted routine
 */
   bet = beta; a = A; b = B; c = C;

   Mjoin( PATL, gegen )( lA, 1, A0, lA, Aseed );
   Mjoin( PATL, gegen )( lB, 1, B0, lB, Bseed );
   Mjoin( PATL, gegen )( lC, 1, C0, lC, Cseed );

#ifdef TREAL
   t0 = time00();
   for( ir = reps; ir; ir-- )
   {
      trusted_symm( SIDE, UPLO, M, N, ALPHA, a, LDA, b, LDB, (SCALAR)(bet),
                    c, LDC );
      a += inca; if( a == stA ) { a = A; }
      b += incb; if( b == stB ) { b = B; }
      c += incc;
      if( c == stC ) { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
   }
   ttrust = time00() - t0;
#else
   if( ROUT == SYMM )
   {
      t0 = time00();
      for( ir = reps; ir; ir-- )
      {
         trusted_symm( SIDE, UPLO, M, N, ALPHA, a, LDA, b, LDB, (SCALAR)(bet),
                       c, LDC );
         a += inca; if( a == stA ) { a = A; }
         b += incb; if( b == stB ) { b = B; }
         c += incc;
         if( c == stC )
         { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
      }
      ttrust = time00() - t0;
   }
   else /* if( ROUT == HEMM ) */
   {
      t0 = time00();
      for( ir = reps; ir; ir-- )
      {
         trusted_hemm( SIDE, UPLO, M, N, ALPHA, a, LDA, b, LDB, (SCALAR)(bet),
                       c, LDC );
         a += inca; if( a == stA ) { a = A; }
         b += incb; if( b == stB ) { b = B; }
         c += incc;
         if( c == stC )
         { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
      }
      ttrust = time00() - t0;
   }
#endif
   if( ttrust > 0.0 ) mftrust = ( reps * flops ) / ( MEGA * ttrust );
   else               mftrust = 0.0;
   ttrust /= reps; *TTRUST0 = ttrust; *MFTRUST0 = mftrust;
   free(A0);
   free(B0);
   free(C0);
   return(0);
}

int syr2kcase
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const int                  MFLOP,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const int                  N,
   const int                  K,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const SCALAR               BETA,
   const int                  LDC,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     flops, ttrust, ttest, mftrust, mftest, t0;
   TYPE                       resid = ATL_rzero;
#ifdef TREAL
   TYPE                       bet,  beta,    nbeta;
#else
   TYPE                       *bet, beta[2], nbeta[2];
#endif
   TYPE                       * a, * stA, *b, * stB, * c, * stC, * A,
                              * A0 = NULL, * B, * B0 = NULL, * C, * C0 = NULL;
   unsigned long              ir, reps;
   int                        inca, incb, incc, lA, lB, lC, mAB, nAB, passed,
                              Aseed, Bseed, Cseed;
   enum ATLAS_TRANS           ta;

   if( ( MEGA * MFLOP <= ( flops = opbl3( ROUT, N, 0, K ) ) ) || ( TEST ) )
   {
      resid = syr2ktst( CACHESIZE, ROUT, TEST, UPLO, TRANS, N, K, ALPHA,
			LDA, LDB, BETA, LDC, EPSILON, TTRUST0, TTEST0,
			MFTRUST0, MFTEST0 );
      if( resid > THRESH ) (void) fprintf( stderr, "   resid=%f\n", resid );
   }
   if( resid < ATL_rzero ) passed = -1;
   else                    passed = ( resid < THRESH );

   if( MEGA * MFLOP <= flops ) return( passed );

   if( TRANS == AtlasNoTrans )
   { ta = TRANS; mAB = N; nAB = K; }
   else
   { ta = ( ROUT == SYR2K ? AtlasTrans : AtlasConjTrans ); mAB = K; nAB = N; }

   inca = LDA  * ( nAB SHIFT );
   incb = LDB  * ( nAB SHIFT );
   incc = LDC  * ( N   SHIFT );

   lA = inca  * ( ( ATL_DivBySize( LCSIZE ) + mAB*nAB - 1 ) / ( mAB * nAB ) );
   lB = incb  * ( ( ATL_DivBySize( LCSIZE ) + mAB*nAB - 1 ) / ( mAB * nAB ) );
   lC = incc  * ( ( ATL_DivBySize( LCSIZE ) + N * N   - 1 ) / ( N   * N   ) );

   A0 = (TYPE *)malloc( ATL_MulBySize( lA ) );
   B0 = (TYPE *)malloc( ATL_MulBySize( lB ) );
   C0 = (TYPE *)malloc( ATL_MulBySize( lC ) );

   if( ( A0 == NULL ) || ( B0 == NULL ) || ( C0 == NULL ) )
   {
      if( A0 ) free( A0 );
      if( B0 ) free( B0 );
      if( C0 ) free( C0 );
      return( -1 );
   }

   A = A0; stA = A0 + ( lA SHIFT );
   B = B0; stB = B0 + ( lB SHIFT );
   C = C0; stC = C0 + ( lC SHIFT );

#ifdef TREAL
   beta   =  BETA;
   nbeta  = -BETA;
#else
   *beta  =    *BETA; beta [1] =  BETA[1];
   *nbeta = -(*BETA); nbeta[1] = -BETA[1];
#endif

   Aseed = mAB * nAB + 513 *  7 + 90;
   Bseed = mAB * nAB + 127 * 50 + 77;
   Cseed = N   * N   + 101 *  2 + 53;

   reps  = ( MEGA * MFLOP ) / flops;
/*
 * Generate the random data and time the trusted routine
 */
   bet = beta; a = A; b = B; c = C;

   Mjoin( PATL, gegen )( lA, 1, A0, lA, Aseed );
   Mjoin( PATL, gegen )( lB, 1, B0, lB, Bseed );
   Mjoin( PATL, gegen )( lC, 1, C0, lC, Cseed );

#ifdef TREAL
   t0 = time00();
   for( ir = reps; ir; ir-- )
   {
      trusted_syr2k( UPLO, ta, N, K, ALPHA, a, LDA, b, LDB, (SCALAR)(bet),
                     c, LDC );
      a += inca; if( a == stA ) { a = A; }
      b += incb; if( b == stB ) { b = B; }
      c += incc;
      if( c == stC ) { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
   }
   ttrust = time00() - t0;
#else
   if( ROUT == SYR2K )
   {
      t0 = time00();
      for( ir = reps; ir; ir-- )
      {
         trusted_syr2k( UPLO, ta, N, K, ALPHA, a, LDA, b, LDB, (SCALAR)(bet),
                        c, LDC );
         a += inca; if( a == stA ) { a = A; }
         b += incb; if( b == stB ) { b = B; }
         c += incc;
         if( c == stC )
         { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
      }
      ttrust = time00() - t0;
   }
   else /* if( ROUT == HER2K ) */
   {
      t0 = time00();
      for( ir = reps; ir; ir-- )
      {
         trusted_her2k( UPLO, ta, N, K, ALPHA, a, LDA, b, LDB, (TYPE)(bet[0]),
                        c, LDC );
         a += inca; if( a == stA ) { a = A; }
         b += incb; if( b == stB ) { b = B; }
         c += incc;
         if( c == stC )
         { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
      }
      ttrust = time00() - t0;
   }
#endif
   if( ttrust > 0.0 ) mftrust = ( reps * flops ) / ( MEGA * ttrust );
   else               mftrust = 0.0;
   ttrust /= reps; *TTRUST0 = ttrust; *MFTRUST0 = mftrust;
   free(A0);
   free(B0);
   free(C0);
   return(0);
}

int syrkcase
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const int                  MFLOP,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const int                  N,
   const int                  K,
   const SCALAR               ALPHA,
   const int                  LDA,
   const SCALAR               BETA,
   const int                  LDC,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     flops, ttrust, ttest, mftrust, mftest, t0;
   TYPE                       resid = ATL_rzero;
#ifdef TREAL
   TYPE                       bet,  beta,    nbeta;
#else
   TYPE                       *bet, beta[2], nbeta[2];
#endif
   TYPE                       * a, * stA, * c, * stC, * A, * A0 = NULL,
                              * C, * C0 = NULL;
   unsigned long              ir, reps;
   int                        inca, incc, lA, lC, mA, nA, passed, Aseed, Cseed;
   enum ATLAS_TRANS           ta;

   if( ( MEGA * MFLOP <= ( flops = opbl3( ROUT, N, 0, K ) ) ) || ( TEST ) )
   {
      resid = syrktst( CACHESIZE, ROUT, TEST, UPLO, TRANS, N, K, ALPHA, LDA,
		       BETA, LDC, EPSILON, TTRUST0, TTEST0, MFTRUST0,
		       MFTEST0 );
      if( resid > THRESH ) (void) fprintf( stderr, "   resid=%f\n", resid );
   }
   if( resid < ATL_rzero ) passed = -1;
   else                    passed = ( resid < THRESH );

   if( MEGA * MFLOP <= flops ) return( passed );

   if( TRANS == AtlasNoTrans )
   { ta = TRANS; mA = N; nA = K; }
   else
   { ta = ( ROUT == SYRK ? AtlasTrans : AtlasConjTrans ); mA = K; nA = N; }

   inca = LDA  * ( nA SHIFT );
   incc = LDC  * ( N  SHIFT );

   lA = inca  * ( ( ATL_DivBySize( LCSIZE ) + mA*nA - 1 ) / ( mA * nA ) );
   lC = incc  * ( ( ATL_DivBySize( LCSIZE ) + N * N - 1 ) / ( N   * N ) );

   A0 = (TYPE *)malloc( ATL_MulBySize( lA ) );
   C0 = (TYPE *)malloc( ATL_MulBySize( lC ) );

   if( ( A0 == NULL ) || ( C0 == NULL ) )
   { if( A0 ) free( A0 ); if( C0 ) free( C0 ); return( -1 ); }

   A = A0; stA = A0 + ( lA SHIFT );
   C = C0; stC = C0 + ( lC SHIFT );

#ifdef TREAL
   beta   =  BETA;
   nbeta  = -BETA;
#else
   *beta  =    *BETA; beta [1] =  BETA[1];
   *nbeta = -(*BETA); nbeta[1] = -BETA[1];
#endif

   Aseed = mA * nA + 513 *  7 + 90;
   Cseed = N  * N  + 101 *  2 + 53;

   reps  = ( MEGA * MFLOP ) / flops;
/*
 * Generate the random data and time the trusted routine
 */
   bet = beta; a = A; c = C;

   Mjoin( PATL, gegen )( lA, 1, A0, lA, Aseed );
   Mjoin( PATL, gegen )( lC, 1, C0, lC, Cseed );

#ifdef TREAL
   t0 = time00();
   for( ir = reps; ir; ir-- )
   {
      trusted_syrk( UPLO, ta, N, K, ALPHA, a, LDA, (SCALAR)(bet), c, LDC );
      a += inca; if( a == stA ) { a = A; }
      c += incc;
      if( c == stC ) { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
   }
   ttrust = time00() - t0;
#else
   if( ROUT == SYRK )
   {
      t0 = time00();
      for( ir = reps; ir; ir-- )
      {
         trusted_syrk( UPLO, ta, N, K, ALPHA, a, LDA, (SCALAR)(bet), c, LDC );
         a += inca; if( a == stA ) { a = A; }
         c += incc;
         if( c == stC )
         { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
      }
      ttrust = time00() - t0;
   }
   else /* if( ROUT == HERK ) */
   {
      t0 = time00();
      for( ir = reps; ir; ir-- )
      {
         trusted_herk( UPLO, ta, N, K, (TYPE)(ALPHA[0]), a, LDA,
                       (TYPE)(bet[0]), c, LDC );
         a += inca; if( a == stA ) { a = A; }
         c += incc;
         if( c == stC )
         { c = C; if( bet == beta ) bet = nbeta; else bet = beta; }
      }
      ttrust = time00() - t0;
   }
#endif
   if( ttrust > 0.0 ) mftrust = ( reps * flops ) / ( MEGA * ttrust );
   else               mftrust = 0.0;
   ttrust /= reps; *TTRUST0 = ttrust; *MFTRUST0 = mftrust;
   free(A0);
   free(C0);
   return(0);
}

int trxmcase
(
   const int                  CACHESIZE,
   const enum LVL3_ROUT       ROUT,
   const int                  TEST,
   const int                  MFLOP,
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_DIAG      DIAG,
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const int                  LDA,
   const int                  LDB,
   const TYPE                 EPSILON,
   double                     * TTRUST0,
   double                     * TTEST0,
   double                     * MFTRUST0,
   double                     * MFTEST0
)
{
   double                     flops, ttrust, ttest, mftrust, mftest, t0;
   TYPE                       resid = ATL_rzero;
   TYPE                       * a, * stA, * b, * stB, * A, * A0 = NULL,
                              * B, * B0 = NULL;
   unsigned long              ir, reps;
   int                        inca, incb, lA, lB, nA, passed, Aseed, Bseed;

   flops = opbl3( ROUT, M, N, ( SIDE == AtlasLeft ? -1 : 1 ) );

   if( ( MEGA * MFLOP <= flops ) || ( TEST ) )
   {
      if( ROUT == TRMM )
      {
         resid = trmmtst( CACHESIZE, ROUT, TEST, SIDE, UPLO, TRANS, DIAG,
			  M, N, ALPHA, LDA, LDB, EPSILON, TTRUST0, TTEST0,
			  MFTRUST0, MFTEST0 );
      }
      else
      {
         resid = trsmtst( CACHESIZE, ROUT, TEST, SIDE, UPLO, TRANS, DIAG,
			  M, N, ALPHA, LDA, LDB, EPSILON, TTRUST0, TTEST0,
			  MFTRUST0, MFTEST0 );
      }
      if( resid > THRESH ) (void) fprintf( stderr, "   resid=%f\n", resid );
   }
   if( resid < ATL_rzero ) passed = -1;
   else                    passed = ( resid < THRESH );

   if( MEGA * MFLOP <= flops ) return( passed );

   if( SIDE == AtlasLeft ) { nA = M; } else { nA = N; }

   inca = LDA * ( nA SHIFT );
   incb = LDB * ( N  SHIFT );

   lA = inca * ( ( ATL_DivBySize( LCSIZE ) + nA*nA - 1 ) / ( nA * nA ) );
   lB = incb * ( ( ATL_DivBySize( LCSIZE ) + M * N - 1 ) / ( M   * N ) );

   A0 = (TYPE *)malloc( ATL_MulBySize( lA ) );
   B0 = (TYPE *)malloc( ATL_MulBySize( lB ) );

   if( ( A0 == NULL ) || ( B0 == NULL ) )
   { if( A0 ) free( A0 ); if( B0 ) free( B0 ); return( -1 ); }

   A = A0; stA = A0 + ( lA SHIFT );
   B = B0; stB = B0 + ( lB SHIFT );

   Aseed = nA * nA + 513 *  7 + 90;
   Bseed = M  * N  + 101 *  2 + 53;

   reps  = ( MEGA * MFLOP ) / flops;
/*
 * Generate the random data and time the trusted routine
 */
   a = A; b = B;

   Mjoin( PATL, gegen )( lA, 1, A0, lA, Aseed );
   Mjoin( PATL, gegen )( lB, 1, B0, lB, Bseed );

   if( ROUT == TRMM )
   {
      t0 = time00();
      for( ir = reps; ir; ir-- )
      {
         trusted_trmm( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, a, LDA, b, LDB );
         a += inca; if( a == stA ) { a = A; }
         b += incb; if( b == stB ) { b = B; }
      }
      ttrust = time00() - t0;
   }
   else /* if( ROUT == TRSM ) */
   {
      do { trddom( UPLO, nA, a, LDA ); a += inca; } while( a != stA ); a = A;

      t0 = time00();
      for( ir = reps; ir; ir-- )
      {
         trusted_trsm( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, a, LDA, b, LDB );
         a += inca; if( a == stA ) { a = A; }
         b += incb; if( b == stB ) { b = B; }
      }
      ttrust = time00() - t0;
   }
   if( ttrust > 0.0 ) mftrust = ( reps * flops ) / ( MEGA * ttrust );
   else               mftrust = 0.0;
   ttrust /= reps; *TTRUST0 = ttrust; *MFTRUST0 = mftrust;
   free(A0);
   free(B0);
   return(0);
}
void PrintUsage(char*, int, char*);
#define ATL_GETFLAGS 1
#include "atlas_genparse.h"

#ifdef TCPLX
   #define NBLAS 9
#else
   #define NBLAS 6
#endif
enum eblas {Bgemm, Bsymm, Bsyrk, Bsyr2k, Btrmm, Btrsm, Bhemm, Bherk, Bher2k};
char *nblas[]={"gemm", "symm", "syrk", "syr2k", "trmm", "trsm", "hemm",
              "herk", "her2k"};

void DoTiming(int nreps, int flszKB, int MFF, int ldaG, int ldbG, int ldcG,
              int ibl, int SD, int UP, int DI, int TA, int TB,
              const SCALAR alpha, const SCALAR beta, int M, int N, int K,
              double *time, double *mflop)
{
   *time = *mflop = 0.0;
   flszKB *= 1024;
   switch(ibl)
   {
      double dtmp;
      int lda, ldb, ldc, ib=0;
   case Bgemm:
      lda = ldaG + (TA == AtlasNoTrans) ? M : K;
      ldb = ldbG + (TB == AtlasNoTrans) ? K : N;
      gemmcase(flszKB, 0, 0, MFF, TA, TB, M, N, K, alpha, lda, M+ldbG, beta,
               M+ldcG, 0.0, time, &dtmp, mflop, &dtmp);
      break;
   case Bhemm:
   case Bsymm:
      ib = (ibl==Bsymm)?SYMM:1;
      lda = ldaG + (SD == AtlasLeft) ? M : N;
      symmcase(flszKB, ib, 0, MFF, SD, UP, M, N, alpha, lda, M+ldbG, beta,
               M+ldcG, 0.0, time, &dtmp, mflop, &dtmp);
      break;
   case Bherk:
   case Bsyrk:
      ib = (ibl==Bsyrk)?SYRK:2;
      lda = ldaG + (TA == AtlasNoTrans) ? N : K;
      syrkcase(flszKB, ib, 0, MFF, UP, TA, N, K, alpha, lda, beta, M+ldcG, 0.0,
               time, &dtmp, mflop, &dtmp);
      break;
   case Bher2k:
   case Bsyr2k:
      ib = (ibl==Bsyr2k)?SYR2K:3;
      lda = ldaG + (TA == AtlasNoTrans) ? N : K;
      ldb = ldbG + (TA == AtlasNoTrans) ? N : K;
      syr2kcase(flszKB, ib, 0, MFF, UP, TA, N, K, alpha, lda, ldb, beta,
                N+ldcG, 0.0, time, &dtmp, mflop, &dtmp);
      break;
   case Btrmm:
   case Btrsm:
      ib = (ibl==Btrsm)?TRSM:TRMM;
      lda = ldaG + (SD == AtlasLeft) ? M : N;
      trxmcase(flszKB, ib, 0, MFF, SD, UP, TA, DI, M, N,
               alpha, lda, M+ldbG, 0.0, time, &dtmp, mflop, &dtmp);
      break;
   }
}

void RunAllTimings(int nreps, int flszKB, int MFF, int ldaG, int ldbG, int ldcG,
                   int *ROUTs, int *SDs, int *UPs, int *DIs, int *TAs, int *TBs,
                   int nalp, TYPE *alps, int nbet, TYPE *bets,
                   int *Ms, int *Ns, int *Ks)
{
   const int nrt = ROUTs[0];
   const char *frm="%6d  %c  %c  %c  %c  %c  %c  %c %7d %7d %7d  %12e %12.1f\n";
   int ir, it=0;

   for (ir=1; ir <= nrt; ir++)
   {
      const int nm = Ms[0];
      int rt = ROUTs[ir], im;
      const int HAS_SD = (rt==Bsymm || rt==Bhemm || rt==Btrmm || rt==Btrsm);
      const int HAS_UP = (rt==Bsyr2k || rt==Bher2k || rt==Bsyrk || rt==Bherk ||
                          rt==Btrmm || rt==Btrsm);
      const int HAS_DI = (rt==Btrmm || rt==Btrsm);
      const int HAS_TA = (rt != Bsymm && rt != Bhemm);
      const int HAS_TB = (rt == Bgemm);
      char cTA='-', cTB='-', cSD='-', cUP='-', cDI='-';

      printf("\n               ********  TIMING FOR %s%-5s ********\n",
             Mstr(PRE), nblas[rt]);
      printf("  TEST SD UP DI TA TB AL BE  "
             "     M       N       K          Time        MFLOP\n");
      printf("====== == == == == == == ==  "
             "======  ======  ======  ============  ===========\n");
      for (im=1; im <= nm; im++)
      {
         const int nn=Ns[0];
         int in;
         for (in=1; in <= nn; in++)
         {
            const int nk = Ks[0];
            const int N = Ns[in], M = (Ms[im]) ? Ms[im] : N;
            int ik;
            for (ik=1; ik <= nk; ik++)
            {
               const int K = (Ks[ik]) ? Ks[ik] : N, ns = SDs[0];
               int is;
               for (is=1; is <= ns; is++)
               {
                  const int SD=SDs[is], nu = UPs[0];
                  int iu;

                  if (HAS_SD)
                     cSD = (SD == AtlasRight) ? 'R' : 'L';
                  for (iu=1; iu <= nu; iu++)
                  {
                     const int UP=UPs[iu];
                     const int nd=DIs[0];
                     int id;

                     if (HAS_UP)
                        cUP = (UP == AtlasLower) ? 'L' : 'U';
                     for (id=1; id <= nd; id++)
                     {
                        const int DI=DIs[id];
                        const int nta=TAs[0];
                        int ita;

                        if (HAS_DI)
                           cDI = (DI == AtlasUnit) ? 'U' : 'N';
                        for (ita=1; ita <= nta; ita++)
                        {
                           const int TA=TAs[ita];
                           const int ntb=TBs[0];
                           int itb;

                           if (HAS_TA)
                           {
                              if (TA == AtlasConjTrans)
                                 cTA = 'C';
                              else
                                 cTA = (TA == AtlasTrans) ? 'T' : 'N';
                           }
                           for (itb=1; itb <= ntb; itb++)
                           {
                              const int TB=TBs[itb];
                              int ia;

                              if (HAS_TB)
                              {
                                 if (TB == AtlasConjTrans)
                                    cTB = 'C';
                                 else
                                    cTB = (TB == AtlasTrans) ? 'T' : 'N';
                              }
                              for (ia=0; ia < nalp; ia++)
                              {
                                 #ifdef TREAL
                                    const TYPE alpha = alps[ia];
                                 #else
                                    const TYPE *alpha = alps+ia+ia;
                                 #endif
                                 int ib;
                                 char cALP;

                                 #ifdef TREAL
                                    if (alpha == 0.0)
                                       cALP = '0';
                                    else if (alpha == 1.0)
                                       cALP = '1';
                                    else if (alpha == -1.0)
                                       cALP = 'N';
                                    else
                                       cALP = 'X';
                                 #else
                                    if (alpha[1] == 0.0)
                                    {
                                       if (*alpha == 0.0)
                                          cALP = '0';
                                       else if (*alpha == 1.0)
                                          cALP = '1';
                                       else if (*alpha == -1.0)
                                          cALP = 'N';
                                       else
                                          cALP = 'R';
                                    }
                                    else
                                       cALP = 'X';
                                 #endif
                                 for (ib=0; ib < nbet; ib++)
                                 {
                                    #ifdef TREAL
                                       const TYPE beta = bets[ia];
                                    #else
                                       const TYPE *beta = bets+ib+ib;
                                    #endif
                                    int r;
                                    char cBET;

                                    #ifdef TREAL
                                       if (beta == 0.0)
                                          cBET = '0';
                                       else if (beta == 1.0)
                                          cBET = '1';
                                       else if (beta == -1.0)
                                          cBET = 'N';
                                       else
                                          cBET = 'X';
                                    #else
                                       if (beta[1] == 0.0)
                                       {
                                          if (*beta == 0.0)
                                             cBET = '0';
                                          else if (*beta == 1.0)
                                             cBET = '1';
                                          else if (*beta == -1.0)
                                             cBET = 'N';
                                          else
                                             cBET = 'R';
                                       }
                                       else
                                          cBET = 'X';
                                    #endif
                                    for (r=0; r < nreps; r++)
                                    {
                                       double time=0.0, mf=0.0;

                                       DoTiming(nreps, flszKB, MFF, ldaG, ldbG,
                                                ldcG, rt, SD, UP, DI, TA, TB,
                                                alpha, beta, M, N, K,
                                                &time, &mf);
                                       printf(frm, ++it, cSD, cUP, cDI, cTA,cTB,
                                              cALP, cBET, M, N, K, time, mf);
                                       fflush(stdout);
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

int *RoutNames2IntList(int nargs, char **args, int i)
{
   int n, *iarr, k;

   if (++i >= nargs)
      PrintUsage(args[0], i, NULL);
   n = atoi(args[i]);
   ATL_assert(n > 0);
   iarr = malloc(sizeof(int)*(n+1));
   ATL_assert(iarr);

   iarr[0] = n;
   for (k=0; k < n; k++)
   {
      int b;
      if (++i >= nargs)
         PrintUsage(args[0], i, NULL);
      for (b=0; b < NBLAS; b++)
      {
         if (!strcmp(args[i], nblas[b]))
         {
            iarr[k+1] = b;
            break;
         }
      }
      if (b == NBLAS)
         PrintUsage(args[0], i, args[i]);
   }
   return(iarr);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -R <#> <rout1> ... <rout#>\n");
   fprintf(stderr, "   -R <rout1/all>\n");
   fprintf(stderr,
   "      routs: gemm, symm, syr2k hemm, her2k, syrk, trmm, trsm, herk\n");
   fprintf(stderr, "   -F <mflop> : force <mflops> of timed computation\n");
   fprintf(stderr, "   -# <#> : repeat each timing # times\n");
   fprintf(stderr, "   -n <#> <N1> ... <N#>\n");
   fprintf(stderr, "   -N <Nstart> <Nend> <Ninc>\n");
   fprintf(stderr, "   -m <#> <M1> ... <M#>\n");
   fprintf(stderr, "   -M <Mstart> <Mend> <Minc>\n");
   fprintf(stderr, "   -k <#> <K1> ... <K#>\n");
   fprintf(stderr, "   -K <Kstart> <Kend> <Kinc>\n");
   fprintf(stderr, "   -G[a,bc] <ldxgap> : ldx = N + <ldxgap> foreach N\n");
   fprintf(stderr, "   -f <flushKB> : flush at least this mem in LRU timers\n");
   fprintf(stderr, "   -S <#> <side1> ... <side#>\n");
   fprintf(stderr, "   -U <#> <up1> ... <up#> : Vals: [u,l,q,g]\n");
   fprintf(stderr, "   -D <#> <dia1> ... <dia#> : Vals: [u,n]\n");
   fprintf(stderr, "   -A <#> <ta1> ... <ta#> : Vals: [n,t,c]\n");
   fprintf(stderr, "   -B <#> <tb1> ... <tb#> : Vals: [n,t,c]\n");
   fprintf(stderr, "   -a <#> <alpha1> ... <alphaN> : real/complex floats\n");
   fprintf(stderr, "   -b <#> <beta1> ... <betaN> : real/complex floats\n");
   exit(ierr ? ierr : -1);
}

void GetFlags(int nargs, char **args, int *nreps, int *flsizeKB, int *mflop,
              int *ldaG, int *ldbG, int *ldcG, int **ROUTs, int **Ns, int **Ms,
              int **Ks, int **UPLOs, int **SDs, int **DIs, int **TAs, int **TBs,
              int *nalphas, TYPE **alphas, int *nbetas, TYPE **betas)
{
   int *NBs=NULL, *ns=NULL, *ms=NULL, *ks=NULL, *ups=NULL, *sds=NULL,
       *dis=NULL, *tas=NULL, *tbs=NULL, *ip;
   int i, j, k, n;

   *ROUTs = NULL;
   *ldbG = *ldcG = *ldaG = 0;
   *flsizeKB = L2SIZE/1024;
   *mflop = 0;
   *nreps = 1;
   *alphas = *betas = NULL;
   for (i=1; i < nargs; i++)
   {
      int WH=0;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'n':         /* -n  or -nb */
         ns = GF_GetIntList(nargs, args, i, 1);
         i += ns[0] + 1;
         break;
      case 'm':         /* -m # <M1> ... <M#> */
         ms = GF_GetIntList(nargs, args, i, 1);
         i += ms[0] + 1;
         break;
      case 'k':         /* -m # <M1> ... <M#> */
         ks = GF_GetIntList(nargs, args, i, 1);
         i += ks[0] + 1;
         break;
      case 'K':         /* -N or -NB */
         WH++;
      case 'N':         /* -N or -NB */
         WH++;
      case 'M':                         /* -M <Mstart> <Mend> <Minc>\n"); */
         if (i+3 >= nargs)
            PrintUsage(args[0], i, NULL);
         ip = GF_IntRange2IntList(atoi(args[i+1]), atoi(args[i+2]),
                                  atoi(args[i+3]));
         if (WH == 2)
            ks = ip;
         else if (WH == 1)
            ns = ip;
         else
            ms = ip;
         i += 3;
         break;
      case 'R':        /* -R # <rout1> ... <routN#>  */
         if (i+1 >= nargs)
            PrintUsage(args[0], i, "out of args to -R");
         if (isdigit(args[i+1][0]))
         {
            *ROUTs = RoutNames2IntList(nargs, args, i);
            i += (*ROUTs)[0] + 1;
         }
         else if (!strcmp(args[++i], "all"))  /* want them all */
            *ROUTs = GF_IntRange2IntList(Bgemm, NBLAS-1, 1);
         else  /* only giving one name */
         {
            for (k=0; k < NBLAS; k++)
               if (!strcmp(args[i], nblas[k]))
                  break;
            if (k == NBLAS)
              PrintUsage(args[0], i,
                         "Unknown blas name to -R: is it lower case?");
            *ROUTs = GF_GetIntList1(k);
         }
         break;
      case '#':                          /* -# <reps> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *nreps = atoi(args[i]);
         break;
      case 'f':                         /* -f <flushKB> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *flsizeKB = atoi(args[i]);
         break;
      case 'F':                         /* -F <mflop> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *mflop = atoi(args[i]);
         break;
      case 'A':                         /* -A <nta> <ta1> ... <taN> */
         WH = 1;
      case 'B':                         /* -B <ntb> <tb1> ... <tbN> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         n = atoi(args[i]);
         ATL_assert(n > 0);
         ip = malloc(sizeof(int)*(n+1));
         ip[0] = n;
         for (k=0; k < n; k++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, NULL);
            switch(args[i][0])
            {
            case 'C':
            case 'c':
               ip[k+1] = AtlasConjTrans;
               break;
            case 't':
            case 'T':
               ip[k+1] = AtlasTrans;
               break;
            default:
               ip[k+1] = AtlasNoTrans;
               break;
            }
         }
         if (WH)
            tas = ip;
         else
            tbs = ip;
         break;
      case 'U':                         /* -U <nup> <u1> ... <uN>;[u,l,q,g] */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         n = atoi(args[i]);
         ATL_assert(n > 0);
         ups = malloc(sizeof(int)*(n+1));
         ups[0] = n;
         for (k=0; k < n; k++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, NULL);
            switch(args[i][0])
            {
            case 'U':
            case 'u':
               ups[k+1] = AtlasUpper;
               break;
            case 'l':
            case 'L':
            default:
               ups[k+1] = AtlasLower;
               break;
            }
         }
         break;
      case 'S':                         /* -S <#> <side1> ... <sideN> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         n = atoi(args[i]);
         ATL_assert(n > 0);
         sds = malloc(sizeof(int)*(n+1));
         sds[0] = n;
         for (k=0; k < n; k++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, NULL);
            switch(args[i][0])
            {
            case 'L':
            case 'l':
               sds[k+1] = AtlasLeft;
               break;
            default:
               sds[k+1] = AtlasRight;
               break;
            }
         }
         break;
      case 'D':                         /* -D <#> <diag1> ... <diagN> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         n = atoi(args[i]);
         ATL_assert(n > 0);
         dis = malloc(sizeof(int)*(n+1));
         dis[0] = n;
         for (k=0; k < n; k++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, NULL);
            switch(args[i][0])
            {
            case 'U':
            case 'u':
               dis[k+1] = AtlasUnit;
               break;
            default:
               dis[k+1] = AtlasNonUnit;
               break;
            }
         }
         break;
      case 'G':                         /* -a <ldagap> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         if (args[i-1][2] == 'b')
            *ldbG = atoi(args[i]);
         else if (args[i-1][2] == 'c')
            *ldcG = atoi(args[i]);
         else
            *ldaG = atoi(args[i]);
         break;
      case 'a':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of flags in -a");
         *nalphas = atoi(args[i]);
         *alphas = malloc(ATL_MulBySize(*nalphas));
         assert(*alphas);
         for (j=0; j < *nalphas SHIFT; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of flags in -a");
            (*alphas)[j] = atof(args[i]);
         }
         break;
      case 'b':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of flags in -b");
         *nbetas  = atoi(args[i]);
         *betas  = malloc(ATL_MulBySize(*nbetas ));
         assert(*betas );
         for (j=0; j < *nbetas SHIFT; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of flags in -b");
            (*betas)[j] = atof(args[i]);
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
/*
 * Take default values
 */
   if (!(*alphas))
   {
      *nalphas = 1;
      *alphas = malloc(sizeof(double)SHIFT);
      **alphas = 1.0;
      #ifdef TCPLX
         (*alphas)[1] = 0.0;
      #endif
   }
   if (!(*betas))
   {
      *nbetas = 1;
      *betas = malloc(sizeof(double)SHIFT);
      **betas = 1.0;
      #ifdef TCPLX
         (*betas)[1] = 0.0;
      #endif
   }
   if (!(*ROUTs))
      *ROUTs = GF_GetIntList1(Bgemm);
   if (!ns)
      ns = GF_GetIntList1(1000);
   if (!ms)
      ms = GF_GetIntList1(0);
   if (!ks)
      ks = GF_GetIntList1(0);
   if (!ups)
      ups = GF_GetIntList1(AtlasUpper);
   if (!sds)
      sds = GF_GetIntList1(AtlasRight);
   if (!dis)
      dis = GF_GetIntList1(AtlasRight);
   if (!tas)
      tas = GF_GetIntList1(AtlasNoTrans);
   if (!tbs)
      tbs = GF_GetIntList1(AtlasNoTrans);

   *Ns = ns;
   *Ms = ms;
   *Ks = ks;
   *UPLOs = ups;
   *SDs = sds;
   *DIs = dis;
   *TAs = tas;
   *TBs = tbs;
}

int main(int nargs, char **args)
{
   int nreps, flszKB, MF, ldaG, ldbG, ldcG, nalp, nbet;
   int *ROUTs, *Ns, *Ms, *Ks, *UPs, *SDs, *DIs, *TAs, *TBs;
   TYPE *alps, *bets;
   GetFlags(nargs, args, &nreps, &flszKB, &MF, &ldaG, &ldbG, &ldcG,
            &ROUTs, &Ns, &Ms, &Ks, &UPs, &SDs, &DIs, &TAs, &TBs,
            &nalp, &alps, &nbet, &bets);
   RunAllTimings(nreps, flszKB, MF, ldaG, ldbG, ldcG, ROUTs, SDs, UPs, DIs,
                 TAs, TBs, nalp, alps, nbet, bets, Ms, Ns, Ks);
   free(ROUTs);
   free(Ms);
   free(Ns);
   free(Ks);
   free(UPs);
   free(SDs);
   free(DIs);
   free(TAs);
   free(TBs);
   free(alps);
   free(bets);
   return(0);
}
