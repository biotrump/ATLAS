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
#include "atlas_pca.h"
#include "atlas_bcamm.h"

static ATL_UINT GetBits(ATL_UINT n)
{
   int i;
   if (n <= 32)
   {
      if (n < 2)
         return(1);
      if (n <= 8) /* 2 <= n <= 8 */
      {
         if (n <= 4) /* 2 <= n <= 4 */
         {
            if (n == 2)
               return(1);
            return(2);
         }
         else        /* 5 <= n <= 8 */
            return(3);
      }
      else if (n <= 16)  /* 8 < n <= 16 */
         return(4);
      else               /* 17 <= n <= 32 */
         return(5);
   }
   for (i=6; (1<<i) < n; i++);
   return(i);
}
/*
 * Decompose a global index into the tuple (rank,lbn,sbn,sbr) where:
 *   rank: rank along given dimension of process grid
 *   lbn : local block number (offset to start of amm's block for index)
 *   sbn : offset within local block of C-format's mu x nu subblock
 *   sbr : row offset with mu x nu sub-block
 * Since each of these quantities is possibly not a power of two, encoding
 * in this form may lose up to 3 bits of range.  This still gives range
 * of 268,435,455 even if we don't use the sign bit.  We will never have
 * enough memory to exceed this for square matrices, and in practice we never
 * will for non-square.  If called with very highly non-square above this
 * range, just use ATLAS's native code rather than calling these amm routs.
 */
void *ATL_bcIpivInit
(
   int MN,         /* length of ipiv array */
   int *ipiv,      /* getrf's ipiv array */
   int inci,       /* piv stride; If inci<0, pivs applied in reverse order */
   int R,          /* max # of procs in row dimension of pgrid */
   int C,          /* max # of procs in column dimension of pgrid */
   int N,          /* dimension of distributed data */
   int B,          /* size of block for distributed dimension */
   int mu,         /* amm unrolling in distributed dim (mu) */
   int nu          /* amm unrolling in non-dist dim (nu) */
)
{
   ATL_UINT nbSBR, nbRNK, nbSBN, n;
   ATL_bcpiv_t *bp;
   void *vp;
   void **lwrks;
   ATL_INT *lldps;
   int i;

/*
 * Return failure if there is not room to safely encode rank in pivot array
 * If this ever happens, simply change ipiv to a 64-bit int, and copy back
 * to 32-bit ipiv at end of computation during decode step.
 */
   nbSBR = GetBits(mu);
   nbRNK = GetBits(R);
   nbSBN = GetBits(B/mu);
   n = 32 - nbSBR - nbRNK - nbSBN;
   if ((1<<n) < N)
      return(NULL);
/*
 * Get required space, initialize it, and return
 */
   bp = malloc(C*sizeof(ATL_bcpiv_t) + MN*sizeof(int) +
                  R*C*(sizeof(void*)+sizeof(int))+ATL_Cachelen);
   if (!bp)
      return(NULL);
   bp->ipiv = (ATL_UINT*)(bp+C); /* for local copy of ipiv0 */
   vp = bp->ipiv + MN;
   lwrks = ATL_AlignPtr(vp);
   lldps = (ATL_INT*)(lwrks+(R*C));

   for (i=0; i<C; i++)
   {
      bp[i].R = R;
      bp[i].C = C;
      bp[i].B = B;
      bp[i].MU = mu;
      bp[i].NU = nu;
      bp[i].neSB = mu*nu;
      bp[i].neMB = B*B;
      bp[i].ipiv = bp->ipiv;
      bp[i].ipiv0 = ipiv;
      bp[i].inci = inci;
      bp[i].nbSBR = nbSBR;
      bp[i].nbRNK = nbRNK;
      bp[i].nbSBN = nbSBN;
      bp[i].larrs = lwrks + (i*R);
      bp[i].lldps = lldps + (i*R);
   }
   return(bp);
}
