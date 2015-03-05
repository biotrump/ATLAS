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
#include "atlas_amm.h"

#ifndef UAMMID
   #define UAMMID 0
#endif
#define PRID Mjoin(Mjoin(atlas_,PRE),Mjoin(Mjoin(u,UAMMID),amm_))
#include Mstr(Mjoin(PRID,blk.h))
#include Mstr(Mjoin(PRID,cm2am_an.h))
#include Mstr(Mjoin(PRID,cm2am_a1.h))
#include Mstr(Mjoin(PRID,cmat2ablk.h))
#include Mstr(Mjoin(PRID,ablk2cmat.h))
#include Mstr(Mjoin(PRID,kern.h))
#ifndef TCPLX
   #define one ATL_rone
   #define none ATL_rnone
   #define zero ATL_rzero
#endif


int Mjoin(PATL,getrf_amm)
  (const int M, const int N, TYPE *A, const int lda, int *ipiv, int IAMM)
{
   int nmb, nnb, mr, nr, NMB, NNB, nb, ierr, minblks, nmu, nnu, mu, nu, m, j;
   size_t sz, csz, dsz, pansz, blksz;
   cm2am_t a2blk, b2blk;
   ablk2cmat_t blk2c;
   cmat2ablk_t c2blk;
   ammkern_t amm;
   void *vp=NULL;
   TYPE *W, *Ac;
   nb = ATL_UAMM_NBs[IAMM];
   b2blk = ATL_UAMM_B2BLK_an[IAMM];
   a2blk = ATL_UAMM_AT2BLK_a1[IAMM];
   c2blk = ATL_UAMM_C2BLK_a1_b0[IAMM];
   blk2c = ATL_UAMM_BLK2C_a1_b0[IAMM];
   amm = ATL_UAMM_KERN_b1[IAMM];
   mu = ATL_UAMM_MUs[IAMM];
   nu = ATL_UAMM_NUs[IAMM];
   if (nb > M || nb > N)
      return(-10);
   nmb = M/nb;
   mr = M - nmb*nb;
   nnb = N/nb;
   nr = N - nnb*nb;
   nmu = nb / mu;
   nnu = nb / nu;
   NMB = (mr) ? nmb+1 : nmb;
   NNB = (nr) ? nnb+1 : nnb;
   if (nr | mr)
      return(-11);  /* handle this later */
   blksz = nb*nb;
   pansz = NMB*blksz;
   sz = NNB*ATL_MulBySize(pansz);
//   if (sz < ATL_MaxMalloc)
      vp = malloc(sz + ATL_Cachelen);
   if (!vp)
      return(-1);
   W = ATL_AlignPtr(vp);
/*
 * Factor first nb-wide panel
 */
   csz = nb*lda;
   dsz = nb*(lda+1);
   ierr = ATL_getrfC(M, nb, A, lda, ipiv);
   Ac = A + csz;
   minblks = Mmin(nmb, nnb);
   m = M - nb;
   for (j=1; j < minblks; j++)
   {
      int i, ie, k;
      TYPE *d = A;                      /* ptr to diag blk to use wt TRSM */
      TYPE *a = Ac;                     /* pointer to blk in apan to update */
      TYPE *apan = W + j*pansz, *b, *c;
      TYPE *l;
/*
 *    Prior iteration has put factored panel back in original storage.
 *    Therefore original storage left of apan has final answer w/o right piv.
 *    Blk storage (W) has A-storage L below diagonal and B-major storage
 *    above diagonal.  The info in W's diagonal blocks should not be used
 *    (it should be solved from A's storage).
 *    At this point, all panels to right and including the active panel
 *    are undefined (not yet copied) in W, and are untouched in original
 *    storage.
 *
 *    For apan, all blocks must have left updates applied.
 *    For the top block, only the solve is necessary, we will do this
 *    in original matrix, then copy it to B-format for U*L of other blocks.
 *    All other blocks must first be copied to C-format so that L on left
 *    can be applied using amm.
 *    After all updates are applied to blocks above diagonal:
 *    1. they are final U (absent right pivot), so copy back to orig storage
 *    2. From there, overwrite C-format storage wt B-format for use in amm.
 *    The diagonal block is copied back to original storage for panel fact.
 *    Blocks below diagonal are copied back to original storage for panel
 *    fact, and then after panel factorization, must be copied to A-format
 *    W for use in applying L (they are read-only after this point).
 */
/*
 *    Apply j*nb pivots to active panel
 *    ERROR HERE.  Can only apply nb at a time, must intermix pivots wt
 *    trsm, I think.
 */
      ATL_laswp(nb, Ac, lda, 0, j*nb, ipiv, 1);
/*
 *    Make first block U in original storage, then copy it to B-format at
 *    top of apan so we can apply first column panel of L.
 */
      cblas_trsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                 CblasUnit, nb, nb, one, d, lda, a, lda);
      b2blk(nb, nb, none, a, lda, apan);
      d += dsz;
      a += nb;    /* finished top block is in orig storage */
/*
 *    All blocks except first have at least one L update to apply, so we
 *    must copy them to C-format in order to apply L*U.
 *    We will apply the first L panel as a special case, since it must
 *    copy from original to C-format, while later L panels can assume
 *    output already in correct format.
 *    This corresponds to a k=0 peeling of following loop, which applies
 *    all L-panels except the first
 */
      c = apan + blksz;
      l = W + blksz;
      for (i=1; i < nmb; i++)   /* first block already final U */
      {
         TYPE *cn = c+blksz, *ln = l+blksz;

         c2blk(nb, nb, one, a, lda, zero, c);
/*
 *       If first thing we apply is just-completed panel fact, must copy it
 *       to A-format!
 */
         if (j == 1)
            a2blk(nb, nb, one, a-csz, lda, l);

         amm(nmu, nnu, nb, l, apan, c, ln, apan, cn);
/*
 *       If this is last update, for entire panel, must copy panel back to
 *       orig storage for use in the panel factorization
 */
         if (j == 1)
            blk2c(nb, nb, one, c, zero, a, lda);
/*
 *       If i==1, then this update frm col pan 0 is all it gets, despite j>1.
 *       Therefore, write it out while in-cache, and then apply TRSM of L
 *       from diag of pan 1 to create final U while in-cache.
 *       Final U must be copied to B-format for use in L*U;
 */
         else if (i == 1)
         {
            blk2c(nb, nb, one, c, zero, a, lda);
            cblas_trsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                       CblasUnit, nb, nb, one, d, lda, a, lda);
            b2blk(nb, nb, none, a, lda, c);
            d += dsz;   /* go to next diagonal block */
         }
         c = cn; l = ln; a += nb;
      }
/*
 *    Now apply remaining L panels to C-storage active panel
 */
      b = apan + blksz;
      for (k=1; k < j; k++)  /* k counts L panels on left of apan */
      {
         l = W + k*(pansz+blksz) + blksz;
         c = apan+k*blksz+blksz;
         a = Ac + k*nb+nb;
         for (i=k+1; i < nmb; i++)
         {
            TYPE *cn = c+blksz, *ln = l+blksz;
/*
 *          If want to apply last update panel (j-1), must copy it to A-form,
 *          since panel factorization left it in original storage
 */
            if (k == j-1)
               a2blk(nb, nb, one, a-csz, lda, l);
            amm(nmu, nnu, nb, l, b, c, ln, b, cn);
/*
 *          If we just applied the last update panel (j-1), we must put
 *          this computed block from active panel back into original
 *          storage for forthcoming use in panel factorization.
 */
            if (k == j-1)
               blk2c(nb, nb, one, c, zero, a, lda);
/*
 *          First block of all panels needs to only have one more L update blk
 *          applied to have finished all updates.  As long there is another
 *          update panel to the left (k+1) of the one we just applied (k)
 *          [this is guaranteed by else], then we can finish this block's
 *          conversion to U while it is still in the cache from the last amm by
 *          copying it back to original storage (a), updating with the k+1
 *          panel's diag block (TRSM) to produce the final U (w/o some pivs).
 *          This is the block immediately above the diagonal block of act pan.
 */
            else if (i == k+1)
            {
               blk2c(nb, nb, one, c, zero, a, lda);
               cblas_trsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                          CblasUnit, nb, nb, one, d, lda, a, lda);
               b2blk(nb, nb, none, a, lda, c);
               d += dsz;   /* go to next diagonal block */
            }
            c = cn; l = ln; a += nb;
         }
         b += blksz;
      }
/*
 *    Perform panel factorization on Lower portion of active panel (orig store)
 */
      i = j*nb;
      ie = ATL_getrfC(m, nb, Ac+i, lda, ipiv+i);
      if (ie && !ierr)
         ierr = ie+i;
/*
 *    Update pivot info to use global index rather than panel-relative.
 *    Normally in left looking we can't apply the new pivots to left, since this
 *    will mess up L, which is read for all remaining cols.  However, we have
 *    copied L into A-major in prior step, it is OK to apply new pivots to L
 *    that has been copied back to original storage, since the algorithm does
 *    not read it.  The only L read from orig storage are diag blocks, which
 *    always have pivots applied via panel fact.
 */
      for (k=i, ie=i+nb; k < ie; k++)
         ipiv[k] += i;
      ATL_laswp(i, A, lda, i, ie-1, ipiv, 1);
      m -= nb;
      Ac += csz;
   }
   free(vp);
   #if 1
   for (j=N-nb-nb; j >= 0; j-=nb)
      ATL_laswp(nb, A+j*nb, lda, j+nb, N, ipiv, 1);
   #endif
   return(ierr);
}
#ifndef TCPLX
   #undef one
   #undef none
   #undef zero
#endif
