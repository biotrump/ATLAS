/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2015 Md Rakib Hasan
 *
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
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
#include "atlas_bcamm.h"

/* copy a panel from cyclic column major to column major */
void Mjoin(PATL, bcL2G_cpy)
(
   int m,
   int n,
   TYPE* W,
   int ldw,
   TYPE* A,
   int lda,
   int nb,
   int nt
)
{
   const int nbnt = nb * nt;
   int m0 = (m/nb)*nb;
   int mr = m - m0;
   int i;
   TYPE *Wc, *Ac=A;
   if (m <= 0) return;
   Ac += ((m0*nt) SHIFT);
   Wc = W + (m0 SHIFT);
   if (mr)
   {
      Mjoin(PATL, gecopy)(mr, n, Wc, ldw, Ac, lda);
   }
   Ac -= (nbnt SHIFT);
   Wc -= (nb SHIFT);
   for (i=m0; i; i-=nb, Ac-=(nbnt SHIFT), Wc-=(nb SHIFT))
   {
      Mjoin(PATL, gecopy)(nb, n, Wc, ldw, Ac, lda);
   }
}
