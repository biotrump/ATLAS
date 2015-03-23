/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.9.24 -- December 25, 2000
 *
 * Author         : Antoine P. Petitet
 * Contributor(s) : R. Clint Whaley
 * Originally developed at the University of Tennessee,
 * Innovative Computing Laboratory, Knoxville TN, 37996-1301, USA.
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the ATLAS group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---------------------------------------------------------------------
 */
/*
 * Include files
 */
#include "atlas_misc.h"
#include "atlas_rblas3.h"
#include "atlas_kernel3.h"
#include "atlas_level3.h"

#ifdef TREAL
   #define H_ALP (const SCALAR)(*((TYPE *)(ALPHA)))
   #define H_BET (const SCALAR)(*((TYPE *)(BETA)))
#else
   #define H_ALP (const SCALAR)(ALPHA)
   #define H_BET (const SCALAR)(BETA)
#endif
/*
 * Type-less wrappers around the ATLAS matrix-multiply functions.
 */
void Mjoin( PATL, gemmNN_RB )
(
   const int                  M,
   const int                  N,
   const int                  K,
   const void                 * ALPHA,
   const void                 * A,
   const int                  LDA,
   const void                 * B,
   const int                  LDB,
   const void                 * BETA,
   void                       * C,
   const int                  LDC
)
{
   Mjoin(PATL,ammm)(AtlasNoTrans, AtlasNoTrans, M, N, K, H_ALP, A, LDA,
                    (const TYPE *)(B), LDB, H_BET, (TYPE *)(C), LDC);
}

void Mjoin( PATL, gemmNT_RB )
(
   const int                  M,
   const int                  N,
   const int                  K,
   const void                 * ALPHA,
   const void                 * A,
   const int                  LDA,
   const void                 * B,
   const int                  LDB,
   const void                 * BETA,
   void                       * C,
   const int                  LDC
)
{
   Mjoin(PATL,ammm)(AtlasNoTrans, AtlasTrans, M, N, K, H_ALP, A, LDA,
                    (const TYPE *)(B), LDB, H_BET, (TYPE *)(C), LDC);
}

void Mjoin( PATL, gemmTN_RB )
(
   const int                  M,
   const int                  N,
   const int                  K,
   const void                 * ALPHA,
   const void                 * A,
   const int                  LDA,
   const void                 * B,
   const int                  LDB,
   const void                 * BETA,
   void                       * C,
   const int                  LDC
)
{
   Mjoin(PATL,ammm)(AtlasTrans, AtlasNoTrans, M, N, K, H_ALP, A, LDA,
                    (const TYPE *)(B), LDB, H_BET, (TYPE *)(C), LDC);
}

#ifdef TCPLX
void Mjoin( PATL, gemmNC_RB )
(
   const int                  M,
   const int                  N,
   const int                  K,
   const void                 * ALPHA,
   const void                 * A,
   const int                  LDA,
   const void                 * B,
   const int                  LDB,
   const void                 * BETA,
   void                       * C,
   const int                  LDC
)
{
   Mjoin(PATL,ammm)(AtlasNoTrans, AtlasConjTrans, M, N, K,
                    (const SCALAR)(ALPHA), (const TYPE *)(A), LDA,
                    (const TYPE *)(B), LDB, (const SCALAR)(BETA),
                    (TYPE *)(C), LDC );
}

void Mjoin( PATL, gemmCN_RB )
(
   const int                  M,
   const int                  N,
   const int                  K,
   const void                 * ALPHA,
   const void                 * A,
   const int                  LDA,
   const void                 * B,
   const int                  LDB,
   const void                 * BETA,
   void                       * C,
   const int                  LDC
)
{
   Mjoin(PATL,ammm)(AtlasConjTrans, AtlasNoTrans, M, N, K,
                    (const SCALAR)(ALPHA), (const TYPE *)(A), LDA,
                    (const TYPE *)(B), LDB, (const SCALAR)(BETA),
                    (TYPE *)(C), LDC );
}
#endif
