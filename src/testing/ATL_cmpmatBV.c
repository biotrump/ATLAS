/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2015 R. Clint Whaley
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
#include "atlas_misc.h"
#include "atlas_bitvec.h"
/*
 * Returns a bitvec where set bits are errors above tolerance tol.
 * Entries are ordered row-major for ease of printing
 * Combine with ATL_print2dBV for pictoral error report.
 */
int *Mjoin(PATL,cmpmatBV)
   (int verb, double tol, int M, int N, const TYPE *A, int lda,
    const TYPE *B, int ldb)
{
   int i, j, *bv;
   size_t lda2=lda SHIFT, ldb2=ldb SHIFT;

   bv = ATL_NewBV(M*N);
   if (tol < 0.0)
      tol = -tol;
   for (j=0; j < N; j++, A += lda2, B += ldb2)
   {
      for (i=0; i < M; i++)
      {
         #ifdef TCPLX
            const int I = i+i;
            double diff = A[I] - B[I], idiff = A[I+1] - B[I+1];
            if (diff < 0.0)
               diff = -diff;
            if (idiff < 0.0)
               idiff = 0.0;
            if (diff > tol || idiff > tol)
            {
               if (verb > 1)
                  printf("A(%d,%d)=[%e,%e];  expected=[%e,%e]\n", i, j,
                         B[I], B[I+1], A[I], A[I+1]);
               ATL_SetBitBV(bv, i*N+j);
            }
         #else
            double diff = A[i] - B[i];
            if (diff < 0.0)
               diff = -diff;

            if (diff > tol)
            {
               if (verb > 1)
                  printf("A(%d,%d)=%e;  expected=%e\n", i, j, B[i], A[i]);
               ATL_SetBitBV(bv, i*N+j);
            }
         #endif

      }
   }
   return(bv);
}
