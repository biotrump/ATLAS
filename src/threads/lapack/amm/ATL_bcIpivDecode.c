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
#include "atlas_pca.h"
#include "atlas_bcamm.h"
/*
 * Takes ipiv split in (rank,lbn, sbn,sbr) bitpattern, and translates it back
 * to standard getrf global row indices.  See ATL_tipivEncode for more details.
 */
void ATL_bcIpivDecode
(
   ATL_bcpiv_t *bp,/* block-cyclic pivot ptr returned by init */
   int n,          /* number of pivot entries to decode */
   int I           /* ipiv index to start the decoding at */
)
{
   ATL_UINT nprior, k;
   ATL_UINT *ipiv = bp->ipiv + I;
   ATL_CUINT P=bp->R, MU=bp->MU, B=bp->B;
   ATL_CUINT nbSBR=bp->nbSBR, nbRNK=bp->nbRNK, nbSBN=bp->nbSBN;
   ATL_CUINT mskSBR=(1<<nbSBR)-1, mskRNK=(1<<nbRNK)-1, mskSBN=(1<<nbSBN)-1;

   for (k=0; k < n; k++)
   {
      ATL_UINT rank, sbn, sbr, lbn = ipiv[k];
      rank = lbn & mskRNK;
      lbn >>= nbRNK;
      sbr = lbn & mskSBR;
      lbn >>= nbSBR;
      sbn = lbn & mskSBN;
      lbn >>= nbSBN;
/*
 *    Scalapack's indxl2g with srcproc=0
 */
      ipiv[k] = lbn*P*B + rank*B + sbn*MU+sbr;
   }
}
