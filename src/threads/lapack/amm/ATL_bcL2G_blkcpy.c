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

/* copy a panel from block-major to cyclic-column major */
void Mjoin(PATL, bcL2G_blkcpy)
(
   int m,
   int n,
   TYPE *W,
   int nb,
   TYPE *A,
   int lda,
   int nt
)
{
   const int nbnb = nb * nb;
   const int nbnt = nb * nt;
   const int m_ = (m/nb)*nb;
   const int mr = m - m_;
   TYPE *Ac, *Wc;
   int i;
   for (i=0, Wc=W, Ac=A; i<m_; i+=nb, Wc+=(nbnb SHIFT), Ac+=(nbnt SHIFT))
   {
      Mjoin(PATL, gecopy)(nb, n, Wc, nb, Ac, lda);
   }
   if (mr)
   {
      Mjoin(PATL, gecopy)(mr, n, Wc, nb, Ac, lda);
   }
}
