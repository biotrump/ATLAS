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
 * ===========================================================================
 * Takes ipiv with lapack-style global row indices, and translates it into
 * block-cyclic friendly tuple  (rank,lbn,sbn,sbr):
 *   rank: prow owning the row the max was found at
 *   lbn: local block number the row is found in
 *   sbn: sub-block in C-major storage row is found in
 *   sbr: sub-block row offset
 * With this tuple, we can compute indxl2g by:
 *   lbn*P*B + rank*B + sbn*mu+sbr
 * And we can determine this local location in the C-format column by:
 *   start = workspace of rank
 *   li = lbn*mb*nb + sbn*mu*nu+sbr
 * Note P (# of procs), B (nb=mb), mu, nu are all stored in bcipiv_t struct.
 *
 * After encoding, ipiv looks like (low order bits on right):
 *   {lbn,sbn,sbr,rank} --> of len --> {32-x, nbSBN, nbSBR, nbRNK}
 *
 *
 * NOTE ON SAFETY OF STORING TUPLE IN ORIGINAL IPIV:
 * -------------------------------------------------
 * For now we store this back in the original ipiv.  In the future we may
 * need to copy it to an N-length 64-bit array, but it should be safe for
 * now.  Note that lbn, sbn, sbr are all indexing the same range as N,
 * so at most we would lose 2 bits of range from this subpartioning
 * (we lose no bits of range if all parts are powers-of-two).  However,
 * the owner rank bits are extra storage, so this is not safe in theory,
 * but it should be in practice.
 * In practice, roughly 1024 cores is the right order to assume for
 * reasonably short term, and this routine is specifically for LU, where
 * we want R<<C (RxC pgrid).  So, this would suggest 32 or 5 bits as
 * a worst-case range loss.  We gain one bit by making unsigned, which
 * gives us 32-5-2+1 = 26 bits of range.  This yields a max N of 65,536
 * for the reasonably foreseeable future.  On current machines with O(100)
 * cores, more like : 32-3-1= 28 = 26,843,546 which is far larger than
 * we can malloc (we have to copy the matrix for this algorithm).
 * So, for now we will return NULL if out of bits (allowing us to safely use
 * ATLAS recursive algorithm, which should be fine performance-wise for
 * asymptotically large problems).  Running out of range is not likely to be
 * a problem before machine evolution forces another threading rewrite.
 * ===========================================================================
 */
/*
 * NOTES:
 * (1) Assumes that all n pivots are in the same block as I.  In practice:
 *     n == nb (< nb for partial blk at end)
 *     I%nb == 0
 * (2) Only the nb rows of IPIV starting at I are encoded.  The destination
 *     rows must be given global index to avoid double encoding!
 */
void ATL_bcIpivEncode
(
   ATL_bcpiv_t *bp,/* block-cyclic pivot ptr returned by init */
   int n,          /* number of pivot entries to encode */
   int I,          /* ipiv index to start the encoding at */
   int iadj        /* amount to adjust pivot entries by */
)
{
   ATL_UINT nprior, k;
   ATL_UINT *ipiv0 = bp->ipiv0 + I;
   ATL_UINT *ipiv = bp->ipiv + I;
   ATL_CUINT P=bp->R, MU=bp->MU, B=bp->B;
   ATL_CUINT nbSBR=bp->nbSBR, nbRNK=bp->nbRNK, nbSBN=bp->nbSBN;
   ATL_CUINT nb2 = nbSBR+nbRNK, nb3=nb2+nbSBN;

   for (k=0; k < n; k++)
   {
      ATL_CUINT d = ipiv0[k]+iadj;    /* global destination row index */
      ATL_CUINT gbn = d/B;           /* global blk num for d'th elt (dst row) */
      ATL_CUINT lbn =  gbn/P;        /* local block # of dest row */
      ATL_CUINT rank = gbn-(lbn)*P;  /* prow of Ith mat row */
      ATL_CUINT boff = d - gbn*B;    /* row within block we want */
      ATL_CUINT sbn = boff/MU;       /* sub-block number of dest row */
      ATL_CUINT sbr = boff - sbn*MU; /* row indx within sub-block of dest row */
/*
 *    Encoding:{lbn,sbn,sbr,rank} --> of len --> {32-x, nbSBN, nbSBR, nbRNK}
 */
      ipiv[k] = (lbn<<nb3) | (sbn<<nb2) | (sbr<<nbRNK) | rank;
   }
}
