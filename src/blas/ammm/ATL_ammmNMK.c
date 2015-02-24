#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_sum.h))

int Mjoin(PATL,ammmNMK)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   const TYPE *B,
   ATL_CSZT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   amminfo_t mminfo;
   ATL_INT mb, NB, kb, mu, nu, ku, KRUN;
   #if ATL_AMM_MAXKMAJ > 1
      size_t kb0U, KK;
   #else
      #define kb0U kb0
      #define KK K
   #endif
   size_t nmblks, nnblks, nkblks, mbF, nbF, kb0, nmu, nmuF, MBF, NBF;
   size_t incwA, incAk, incBk, incAk0, incBk0, incAm, incBn, incCn;
   size_t i, j, k, szA, szB, szC, nnu0, nnuF;
   const TYPE *B0 = B;
   TYPE alpA=ATL_rone, alpB=ATL_rone, alpC=ATL_rone;
   TYPE *wA, *wB, *wC, *wA0, *wB0;
   void *vp;
   cm2am_t a2blk, b2blk;
   ablk2cmat_t blk2c;
   ammkern_t amm_b1, amm_b0;
   int appAl;

   appAl = Mjoin(PATL,GetAmmmInfo)(&mminfo, TA, TB, M, N, K, alpha, beta);
   if (!appAl)
      alpA = alpha;
   else if (appAl == 1)
      alpB = alpha;
   else
      alpC = alpha;
   mb = mminfo.mb;
   NB = mminfo.nb;
   kb = mminfo.kb;
   KRUN = ATL_AMMFLG_KRUNTIME(mminfo.flag);
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   incwA = mb*kb;
   nmu = mb / mu;
   amm_b1 = mminfo.amm_b1;
   blk2c = mminfo.Cblk2cm;
   a2blk = mminfo.a2blk;
   b2blk = mminfo.b2blk;
/*
 * Handle N differently than other dims: since it is outer loop, don't want
 * to peel, so just vary nb inside the main loop.  nnblks therefore includes
 * final block, unlike for M or K.
 */
   if (N >= NB+nu+nu)
   {
      nnblks = N/NB;
      nbF = N - nnblks * NB;
      if (nbF < nu+nu)
         nbF += NB;
      else
         nnblks++;
   }
   else
   {
      nnblks = 1;
      nbF = N;
   }
   nnu0 = NB / nu;
   nnuF = (nbF+nu-1)/nu;
   NBF = nnuF * nu;
/*
 * For M, we must peel the final block to handle any cleanup (can't peel 1st
 * block or we mess up alignment!), so this block is not included in the
 * block count
 */
   if (M >= mb+mu+mu)  /* more than just last block */
   {
      nmblks = M/mb;
      mbF = M - nmblks * mb;
      if (mbF < mu+mu)  /* steal block from main iteration, not enough here! */
      {
         nmblks--;
         mbF += mb;
      }
   }
   else /* put everything in final block */
   {
      mbF = M;
      nmblks = 0;
   }
   nmuF = (mbF+mu-1) / mu;
   MBF = nmuF * mu;
/*
 * For K, we peel the first iteration to set BETA=0, so the nkblks does not
 * include the peeled block
 */
   if (K >= kb)
   {
      nkblks = K/kb;
      kb0 = K - nkblks * kb;
      if (!kb0)
      {
         kb0 = kb;
         nkblks--;
      }
      else
      {
         if (kb0 < 4)
         {
            kb0 += kb;
            nkblks--;
         }
      }
   }
   else /* K < nb */
   {
      kb0 = K;
      nkblks = 0;
   }
   #if ATL_AMM_MAXKMAJ > 1
      if (ATL_AMMFLG_KMAJOR(mminfo.flag))
      {
         KK = ((K+ku-1)/ku)*ku;
         kb0U = ((kb0+ku-1)/ku)*ku;
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag))
            amm_b0 = mminfo.amm_b0;
         else
            amm_b0 = (mminfo.kb == kb0U) ?  mminfo.amm_b0 : mminfo.amm_k1_b0;
      }
      else
      {
         KK = K;
         kb0U = kb0;
         amm_b0 = (kb0 == mminfo.kb ||
                   (KRUN && kb0 >= mminfo.kbmin && (kb0/ku)*ku == kb0)) ?
                  mminfo.amm_b0 : mminfo.amm_k1_b0;
      }
   #else
      amm_b0 = (kb0 == mminfo.kb ||
                (KRUN && kb0 >= mminfo.kbmin && (kb0/ku)*ku == kb0)) ?
               mminfo.amm_b0 : mminfo.amm_k1_b0;
   #endif
   szA = (nmblks*mb+MBF)*KK; /* wrkspc for all of A wt M rounded up to MU*/
   j = Mmax(NB, NBF);
   i = Mmax(MBF, mb);
   szC = i*j;
   szB = KK*j;                     /* workspace for panel of B */

   k = ATL_MulBySize(szA+szB+szC+mu*nu*ku) + 3*ATL_Cachelen;
   if (k > ATL_MaxMalloc)
      return(2);
   vp = malloc(k);
   if (!vp)
      return(1);
   wB0 = wB = ATL_AlignPtr(vp);
   wA = wB + szB;
   wA0 = wA = ATL_AlignPtr(wA);
   wC = wA + szA;
   wC = ATL_AlignPtr(wC);

   if (TA == AtlasNoTrans)
   {
      incAm = mb;
      incAk0 = kb0*lda;
      incAk = kb*lda;
   }
   else
   {
      incAm = mb*lda;
      incAk0 = kb0;
      incAk = kb;
   }
   if (TB == AtlasNoTrans)
   {
      incBk0 = kb0;
      incBk = kb;
      incBn = NB*ldb;
   }
   else
   {
      incBk0 = kb0*ldb;
      incBk = kb*ldb;
      incBn = NB;
   }
   incCn = ldc*NB;

   for (j=0; j < nnblks; j++)
   {
      size_t nb, nbsz, incwB, incwB0, nnu;
      const TYPE *Bn = B+incBn;
      TYPE *Cn = C + incCn;
      if (j != nnblks-1)
      {
         nbsz = nb = NB;
         nnu = nnu0;
      }
      else
      {
         nb = nbF;
         nbsz = NBF;
         nnu = nnuF;
      }
      incwB = nbsz*kb;
      incwB0 = nbsz*kb0U;
/*
 *    Do all M-blocks except final one, which may be of differing size & partial
 */
      for (i=0; i < nmblks; i++)
      {
         TYPE *wAn, *wBn;
         const TYPE *An = A+incAm;
/*
 *       Peel first K it to handle K-cleanup and set BETA=0
 */
         if (!j)
            a2blk(kb0, mb, alpA, A, lda, wA);
         if (!i)
            b2blk(kb0, nb, alpB, B, ldb, wB);
         wAn = wA+mb*kb0U;
         wBn = (nkblks) ? wB+incwB0 : wB;
         amm_b0(nmu, nnu, kb0U, wA, wB, wC, nkblks?wAn:wA, wBn, wC);
         wA = wAn;
         wB = wBn;
         A += incAk0;
         B += incBk0;
/*
 *       If first K-block not the only K-block
 */
         if (nkblks)
         {
            for (k=nkblks-1; k; k--)
            {
               wAn = wA+incwA;
               wBn = wB+incwB;
               if (!j)
                  a2blk(kb, mb, alpA, A, lda, wA);
               if (!i)
                  b2blk(kb, nb, alpB, B, ldb, wB);
               amm_b1(nmu, nnu, kb, wA, wB, wC, wAn, wBn, wC);
               wA = wAn;
               wB = wBn;
               A += incAk;
               B += incBk;
            }
/*
 *          Last K-block peeled to change prefetch pattern
 */
            wAn = wA+incwA;
            if (!j)
               a2blk(kb, mb, alpA, A, lda, wA);
            if (!i)
               b2blk(kb, nb, alpB, B, ldb, wB);
            amm_b1(nmu, nnu, kb, wA, wB, wC, wAn, wB0, wC);
            wA = wAn;
            wB = wB0;
         }
         blk2c(mb, nb, alpC, wC, beta, C, ldc);
         A = An;
         B = B0;
         C += mb;
      }
/*
 *    Do the final peeled M-block, which is of non-constant size mbF
 */
      {
         TYPE *wAn, *wBn;
/*
 *       Peel first K it to handle K-cleanup and set BETA=0
 */
         if (!j)
            a2blk(kb0, mbF, alpA, A, lda, wA);
         if (!i)
            b2blk(kb0, nb, alpB, B, ldb, wB);
         wAn = wA+MBF*kb0U;
         wBn = (nkblks) ? wB+incwB0 : wB;
         amm_b0(nmuF, nnu, kb0U, wA, wB, wC, nkblks?wAn:wA, wBn, wC);
         wA = wAn;
         wB = wBn;
         A += incAk0;
         B += incBk0;
/*
 *       If first K-block not the only K-block
 */
         if (nkblks)
         {
            for (k=nkblks-1; k; k--)
            {
               wAn = wA+MBF*kb;
               wBn = wB+incwB;
               if (!j)
                  a2blk(kb, mbF, alpA, A, lda, wA);
               if (!i)
                  b2blk(kb, nb, alpB, B, ldb, wB);
               amm_b1(nmuF, nnu, kb, wA, wB, wC, wAn, wBn, wC);
               wA += MBF*kb;
               wB = wBn;
               A += incAk;
               B += incBk;
            }
/*
 *          Last K-block peeled to change prefetch pattern
 */
            if (!j)
               a2blk(kb, mbF, alpA, A, lda, wA);
            if (!i)
               b2blk(kb, nb, alpB, B, ldb, wB);
            amm_b1(nmuF, nnu, kb, wA, wB, wC, wA0, wB0, wC);
            wA = wA0;
            wB = wB0;
         }
         blk2c(mbF, nb, alpC, wC, beta, C, ldc);
      }  /* end M-peel */
      wA = wA0;
      B = B0 = Bn;
      C = Cn;
   }
   free(vp);
   return(0);
}
