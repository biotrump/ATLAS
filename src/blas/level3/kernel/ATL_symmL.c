/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 *                    (C) Copyright 1997 R. Clint Whaley
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
#include "atlas_kern3.h"

#ifdef TREAL
   void Mjoin(Mjoin(Mjoin(PATL,sycopy),UploNM),_a1)
      (const int N, const SCALAR alpha0, const TYPE *A, const int lda, TYPE *C);
   void Mjoin(Mjoin(Mjoin(PATL,sycopy),UploNM),_aX)
      (const int N, const SCALAR alpha0, const TYPE *A, const int lda, TYPE *C);
#else
   void Mjoin(Mjoin(PATL,sycopy),UploNM)
      (const int N, const TYPE *A, const int lda, TYPE *C);
#endif

void Mjoin(Mjoin(PATL,symmL),UploNM)
   (const int M, const int N, const void *valpha, const void *A, const int lda,
    const void *B, const int ldb, const void *vbeta, void *C, const int ldc)
{
   #ifdef TREAL
      const SCALAR alpha=*( (const SCALAR *)valpha );
      const SCALAR beta =*( (const SCALAR *)vbeta  );
      const SCALAR one=1.0;
   #else
      #define alpha valpha
      #define beta vbeta
   #endif
   TYPE *a;
   void *va;

   if (N > SYMM_Xover)
   {
      va = malloc(ATL_Cachelen + (ATL_MulBySize(M)*M));
      ATL_assert(va);
      a = ATL_AlignPtr(va);
      #ifdef TREAL
         if ( SCALAR_IS_ONE(alpha) )
            Mjoin(Mjoin(Mjoin(PATL,sycopy),UploNM),_a1)(M, alpha, A, lda, a);
         else Mjoin(Mjoin(Mjoin(PATL,sycopy),UploNM),_aX)(M, alpha, A, lda, a);
         ATL_gemm(AtlasTrans, AtlasNoTrans, M, N, M, one, a, M, B, ldb,
                  beta, C, ldc);
      #else
         Mjoin(Mjoin(PATL,sycopy),UploNM)(M, A, lda, a);
         ATL_gemm(AtlasTrans, AtlasNoTrans, M, N, M, valpha, a, M, B, ldb,
                  vbeta, C, ldc);
      #endif
      free(va);
   }
   else Mjoin(PATL,refsymm)(AtlasLeft, Uplo_, M, N, alpha, A, lda, B, ldb,
                            beta, C, ldc);
}
