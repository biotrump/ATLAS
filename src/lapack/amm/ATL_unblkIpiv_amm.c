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
#include "atlas_lapack.h"
/*
 * Takes ipiv split in (mbn,sbn,sbr) bitpattern, and translates it back
 * to standar getrf row coordinates.  See ATL_blkIPiv_amm for more detail
 */
void ATL_unblkIpiv_amm
(
   int n,      /* number of pivot entries to unblock */
   int *ipiv,  /* pivot array to encode wt block bit patterns */
   const int nb,
   const int mu,
   const int bb,
   const int sb,
   const int rb
)
{
   unsigned int k;
   const unsigned int bmsk=(1<<bb)-1, smsk=(1<<sb)-1, bb_sb=bb+sb;
   for (k=0; k < n; k++)
   {
      const unsigned int i = ipiv[k];
      unsigned int mbn, sbn, sbr;
      mbn = i & bmsk;
      sbn = (i>>bb) & smsk;
      sbr = (i>>bb_sb);
      ipiv[k] = sbr + sbn*mu + mbn*nb;
//printf("decode: (%d,%d,%d) = %d\n", mbn, sbn, sbr, ipiv[k]);
   }
}

