/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2015, 2013, 2012 R. Clint Whaley
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
#include "atlas_gnuvec.h"
#define ATL_JKMDEF 1
#include "atlas_mmtesttime.h"

#define NVECS 4
static enum VECTYPE {VTAVXZ=0, VTAVX=1, VTSSE=2, VTGV=3, VTSC=4} VECi=VTSC;
static int VLEN[5] = {8, 4, 2, 2, 1};  /* assume double, fix later if nec */
static char *VECs[5] = {"avxz", "avx", "sse", "gvec", "scalar"};
static int TSIZE=8;
static int IMVS=3;     /* move ptrs in timing encoded in last 3 bits: CBA */
#define KRUNMUL 1.02   /* KRUNTIME speedup increase over K-compile time */

static int Mylcm(const int M, const int N)
/*
 * Returns least common multiple (LCM) of two positive integers M & N by
 * computing greatest common divisor (GCD) and using the property that
 * M*N = GCD*LCM.
 */
{
   register int tmp, max, min, gcd=0;

   if (M != N)
   {
      if (M > N) { max = M; min = N; }
      else { max = N; min = M; }
      if (min > 0)  /* undefined for negative numbers */
      {
         do  /* while (min) */
         {
            if ( !(min & 1) ) /* min is even */
            {
               if ( !(max & 1) ) /* max is also even */
               {
                  do
                  {
                     min >>= 1;
                     max >>= 1;
                     gcd++;
                     if (min & 1) goto MinIsOdd;
                  }
                  while ( !(max & 1) );
               }
               do min >>=1 ; while ( !(min & 1) );
            }
/*
 *          Once min is odd, halve max until it too is odd.  Then, use
 *          property that gcd(max, min) = gcd(max, (max-min)/2)
 *          for odd max & min
 */
MinIsOdd:
            if (min != 1)
            {
               do  /* while (max >= min */
               {
                  max -= (max & 1) ? min : 0;
                  max >>= 1;
               }
               while (max >= min);
            }
            else return( (M*N) / (1<<gcd) );
            tmp = max;
            max = min;
            min = tmp;
         }
         while(tmp);
      }
      return( (M*N) / (max<<gcd) );
   }
   else return(M);
}

double TimeMMKernel_KB
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int FORCETIME,               /* 1: ignore any prior output file */
   ATL_mmnode_t *mmp,           /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
/*
 * If kernel has property KRUNTIME, try timing it with compile- and run-time K,
 * and if compile-time is more than 2% faster, turn off KRUNTIME
 */
{
   double mf;
   if (mmp->kbmax && kb > mmp->kbmax)
      return(0.0);
   if (mmp->kbmin && kb < mmp->kbmin)
      return(0.0);
   mf = TimeMMKernel(verb, FORCETIME, mmp, pre, mb, nb, kb, 0, 0, 0,
                     beta, mflop, cflush);
   if (FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
   {
      double mfC;
      mmp->flag &= ~(1<<MMF_KRUNTIME);
      mfC = TimeMMKernel(verb, FORCETIME, mmp, pre, mb, nb, kb, kb, kb, mb,
                         beta, mflop, cflush);
      if (mfC <= 1.02*mf)
         mmp->flag |= (1<<MMF_KRUNTIME);
      else
      {
         if (verb)
            printf("      Forcing K compile-time, mfC=%.2f, mfR=%.2f\n",
                   mfC, mf);
         mf = mfC;
      }
   }
   return(mf);
}
/*
 * Finds best blocking factors for kernel mmp trying all legal values
 * between [b0, bN]
 */
ATL_mmnode_t *BestBlocking_BFI
(
   int verb,
   char pre,
   ATL_mmnode_t *mmp,
   int b0,
   int bN,
   int minInc,  /* minimum increment to use */
   int FORCE
)
/*
 * Times all legal block factors in all dims between [b0,bN].
 * RETURNS: ptr to best performing kernel, NULL if no legal block factors
 */
{
   ATL_mmnode_t *mp;
   int mbB=0, nbB=0, kbB=0;
   int mbS=0, nbS=0, kbS=0;
   int mu = mmp->mu, nu = mmp->nu, ku = mmp->ku;
   int k0, kn, m0, mn, n0, nn, m, n, k;
   double mfB=0.0, mfS=0.0;

   if (!mmp)
      return(NULL);
   if (minInc > mu)
      mu = ((minInc+mu-1)/mu)*mu;
   if (minInc > nu)
      nu = ((minInc+nu-1)/nu)*nu;
   if (minInc > ku)
      ku = ((minInc+ku-1)/ku)*ku;
   m0 = ((b0+mu-1)/mu)*mu;
   n0 = ((b0+nu-1)/nu)*nu;
   k0 = ((b0+ku-1)/ku)*ku;
   mn = ((bN+mu-1)/mu)*mu;
   nn = ((bN+nu-1)/nu)*nu;
   kn = ((bN+ku-1)/ku)*ku;
   mp = CloneMMNode(mmp);
   if (mp->kbmax && mp->kbmax < kn)
      kn = mp->kbmax;
   if (mp->kbmin && mp->kbmin > k0)
      k0 = mp->kbmin;


   printf("SEARCH BLKING [%d - %d] for %d.%s:\n\n", b0, bN, mp->ID, mp->rout);
   printf("  MB    NB    KB        MFLOP    mbB  nbB  kbB      mflopB\n");
   printf("====  ====  ====  ===========   ==== ==== ==== ===========\n");
   for (m=m0; m <= mn; m += mu)
   {
      for (n=m0; n <= nn; n += nu)
      {
         for (k=k0; k <= kn; k += ku)
         {
            double mf;
            mf = TimeMMKernel(verb, FORCE, mp, pre, m, n, k, k, k, m, 1, 0, -1);
            printf("%4d %5d %5d %11.1f %4d %4d %4d %11.1f\n",
                   m, n, k, mf, mbB, nbB, kbB, mfB);
            if (mf > mfB)
            {
               mfB = mf;
               mbB = m;
               nbB = n;
               kbB = k;
            }
            if (m == n && m == k)
            {
               if (mf > mfS)
               {
                  mfS = mf;
                  mbS = m;
                  nbS = n;
                  kbS = k;
               }
            }
         }
      }
   }
   if (mfB == 0)
   {
      printf("NO KERNEL POSSIBLE FOR RANGE=[%d,%d]\n", b0, bN);
      KillMMNode(mp);
      return(NULL);
   }
   mp->mbB = mbB;
   mp->nbB = nbB;
   mp->kbB = kbB;
   mp->mflop[0] = mfB;
   printf("FOR %d.'%s': MB=%d, NB=%d, KB=%d, MFLOPS=%.1f\n",
          mp->ID, mp->rout, mbB, nbB, kbB, mfB);
   k = MMKernelFailsTest(pre, mbB, nbB, kbB, 0, mp);
   if (!k)
      k = MMKernelFailsTest(pre, mbB, nbB, kbB, 1, mp);
   if (!k)
      k = MMKernelFailsTest(pre, mbB, nbB, kbB, -1, mp);
   if (k)
   {
      printf("KERNEL FAILS TESTER FOR [M,N,K]B=%d,%d,%d\n", mbB, nbB, kbB);
      exit(k);
   }
   if (mbS == 0)
      mp->next = NULL;
   else
   {
      k = MMKernelFailsTest(pre, mbS, nbS, kbS, 0, mp);
      if (!k)
         k = MMKernelFailsTest(pre, mbS, nbS, kbS, 1, mp);
      if (!k)
         k = MMKernelFailsTest(pre, mbS, nbS, kbS, -1, mp);
      if (k)
         mp->next = NULL;
      else
      {
         mp->next = CloneMMNode(mp);
         mp->next->mbB = mbS;
         mp->next->nbB = nbS;
         mp->next->kbB = kbS;
         mp->next->mflop[0] = mfS;
      }
   }
   WriteMMFileWithPath(pre, "res", "AMMEXBLKS.sum", mp);
   return(mp);
}

ATL_mmnode_t *TimeExtraBlockings(char pre, int verb)
{
   ATL_mmnode_t *eb;
   eb = ReadMMFileWithPath(pre, "res", "AMMEXBLKS.sum");
   if (!eb)
      return(eb);
   if (eb->mflop[0] < 0)
   {
      ATL_mmnode_t *mp;
      printf("EXTRA BLOCKING FACTOR TIMINGS:\n\n");
      if (verb)
      {
         printf("  MB    NB    KB        MFLOP\n");
         printf("====  ====  ====  ===========\n");
      }
      for (mp=eb; mp; mp = mp->next)
      {
         mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mp->mbB, mp->nbB,
                                     mp->kbB, 0, 0, mp->mbB, 1, 0, -1);
         if (verb)
            printf("%4d %5d %5d %11.1f\n",
                   mp->mbB, mp->nbB, mp->kbB, mp->mflop[0]);
      }
      WriteMMFileWithPath(pre, "res", "AMMEXBLKS.sum", eb);
   }
   return(eb);
}

ATL_mmnode_t *BestForThisNB
(
   int verb,
   char pre,
   ATL_mmnode_t *mmb,
   int nb,
   int pnb,  /* previous nb */
   int nnb,  /* next nb */
   int FORCE
)
/*
 * Times all kernels in mmb
 * RETURNS: ptr to best performing kernel, empty gen node if no user case wrks
 */
{
   ATL_mmnode_t *mmp, *mmB=NULL;
   double mf, mf0, mfB=0.0;

   printf("SCOPING FOR BEST PERFORMING KERNEL FOR NB=%d\n", nb);
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      int kb;
/*
 *    Choose kb, if forced only kb will do, so skip if kernel can't do it
 */
      if (FORCE || nb <= 16)
         kb = nb;
/*
 *    If this kernel can't do the exact block factor, allow leeway
 */
      else
      {
         int u;
         u = Mylcm(mmp->mu, mmp->nu);
         u = Mylcm(u, mmp->ku);
         kb = (nb/u)*u;
         if (kb != nb)
         {
            int kbB;
            kbB = ((nb+u-1)/u)*u;
            if (!kb || kbB-nb < nb-kb && kbB <= 4)
               kb = kbB;
         }
      }
      if ((mmp->kbmin && kb < mmp->kbmin) ||
          (mmp->kbmax && kb > mmp->kbmax) ||
          ((kb/mmp->mu)*mmp->mu != kb) || ((kb/mmp->nu)*mmp->nu != kb) ||
          ((kb/mmp->ku)*mmp->ku != kb) || (kb == pnb) || (kb == nnb))
      {

         printf("   %d. %s: SKIPPED, bad NB\n", mmp->ID, mmp->rout);
         continue;
      }
      mf0 = TimeMMKernel(verb, 0, mmp, pre, kb, kb, kb, kb, kb, kb, 1, 0, -1);
/*
 *    Give bonus to K-runtime variable over K-compile time; K-runtime kernels
 *    can be used for some K-cleanup, and they can be used for any required KB
 *    as well as being typically much smaller instruction load, so they are
 *    strongly preferred
 */
      mf = FLAG_IS_SET(mmp->flag, MMF_KRUNTIME) ? mf0*KRUNMUL : mf0;
      if (mf > mfB)
      {
         mfB = mf;
         mmB = mmp;
         mmB->mbB = mmB->nbB = mmB->kbB = kb;
      }
      printf("   %d. %s: kb=%d, MFLOP=%.2f\n", mmp->ID, mmp->rout, kb, mf0);
   }
   if (!mmB)
   {
      printf("NO KERNEL POSSIBLE FOR NB=%d\n", nb);
      mmB = GetMMNode();
      mmB->mbB = mmB->nbB = mmB->kbB = nb;
   }
   else
   {
      int i, kb = mmB->kbB;
      i = MMKernelFailsTest(pre, kb, kb, kb, 0, mmB);
      if (!i)
         i = MMKernelFailsTest(pre, kb, kb, kb, 1, mmB);
      if (!i)
         i = MMKernelFailsTest(pre, kb, kb, kb, -1, mmB);
      if (i)
      {
         printf("BEST KERNEL FAILS TESTER FOR NB=%d\n", kb);
         exit(i);
      }
      if (FLAG_IS_SET(mmB->flag, MMF_KRUNTIME))
         mfB /= KRUNMUL;
      printf("BEST KERNEL FOUND FOR NB=%d: ID#%d '%s' %.2f MFLOPS\n\n",
             nb, mmB->ID, mmB->rout, mfB);
      mmB = CloneMMNode(mmB);
      mmB->mflop[0] = mfB;
      mmB->next = NULL;
   }
   return(mmB);
}

int DeleteBadBigNBs(ATL_mmnode_t *mmb, int *nbs)
{
   ATL_mmnode_t *best=NULL, *mmp;
   double mfB=0.0;
   int n=0;
/*
 * Find the best-performing kernel
 */
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      double mf;
      mf = mmp->mflop[0];
      if (mf > mfB)
      {
         mfB = mf;
         best = mmp;
      }
   }
/*
 * Delete all NBs larger than best
 */
   while (best->next)
   {
      best->next = KillMMNode(best->next);
      n++;
   }
   if (n)
   {
      int N = *nbs;
      N = (N >= 0) ? N : -N;
      printf("Deleted %d large, slow kernels starting at NB=%d\n",
             n, nbs[N-n+1]);
   }
   return(n);
}

ATL_mmnode_t *FindBestForEachNB(int verb, char pre, ATL_mmnode_t *mmb, int *nbs)
{
   int i, n, FORCE=0;
   ATL_mmnode_t *best, *bp;
/*
 * If # of nbs is negative, then each nb is required and that exact size
 * will be used, or no NB of that size if no kernel works.  The normal behavior
 * is the exact size of forced for all nb <= 16, and inexact for larger
 */
   n = nbs[0];
   if (n < 0)  /* negative # of nbs says force exact NB or nothing */
   {
      n = -n;
      FORCE = 1;
   }
   bp = best = BestForThisNB(verb, pre, mmb, nbs[1], 0, (n == 1)?nbs[1]:0,
                             FORCE);
   for (i=2; i <= n; i++)
   {
      int pnb = nbs[i-1], nnb = (i < n) ? nbs[i+1]:0;
      bp->next = BestForThisNB(verb, pre, mmb, nbs[i], pnb, nnb, FORCE);
      bp = bp->next;
   }
   if (!FORCE)
      i = DeleteBadBigNBs(best, nbs);
   return(best);
}

static INLINE void ApplyMoves2Flag
(
   ATL_mmnode_t *mmp,  /* kernel to set MMF_MV[A,B,C] flag bits */
   int mvBits          /* last 3 bits: MOVE_[CBA] */
)
{
   int flag = mmp->flag & (~MMF_MVSET);         /* zero existing move bits */
   mmp->flag = flag | ((mvBits & 7)<<MMF_MVA); /* put new move bits in */
}
static void ApplyMoves2Flags
(
   ATL_mmnode_t *mmb,  /* kernel to set MMF_MV[A,B,C] flag bits */
   int mvBits          /* last 3 bits: MOVE_[CBA] */
)
{
   const unsigned int mvMSK = ~MMF_MVSET, mvSET = (mvBits&7)<<MMF_MVA;
   ATL_mmnode_t *mmp;
   for (mmp=mmb; mmp; mmp = mmp->next)
      mmp->flag = ((mmp->flag) & mvMSK) | mvSET;
}

char *GenString(char pre, int lat, int nb, int mu, int nu, int ku,
                int kmaj, char *rt)
{
   char *frm=
      "make gen_amm_%s lat=%d mu=%d nu=%d ku=%d kb=%d vlen=%d rt=%s kmaj=%d";
   char *ln;
   int l;
   l = strlen(frm) + strlen(VECs[VECi]) + strlen(rt) + 8;
   ln = malloc(l*sizeof(char));
   assert(ln);
   if (kmaj > 1)
      sprintf(ln, frm, "scalar", lat, mu, nu, ku, nb, 1, rt, kmaj);
   else if (VECi != VTGV)
      sprintf(ln, frm, VECs[VECi], lat, mu*VLEN[VECi], nu, ku, nb,
              VLEN[VECi], rt, 0);
   else
   {
      int sz = (pre == 's') ? 4 : 8;
      sprintf(ln, frm, VECs[VECi], lat, mu, nu, ku, nb, VLEN[VECi]*sz, rt, 0);
   }
   return(ln);
}

void FillInGenStrings(char pre, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mp;
   for (mp=mmb; mp; mp = mp->next)
   {
      if (mp->ID == 0 && !mp->genstr)
      {
         const int vlen=mp->vlen, mu = (vlen) ? mp->mu/vlen : mp->mu, vl0=VECi;
         if (!vlen || vlen == 1)
            VECi = VTSC;
         else if (vlen != VLEN[VECi])
         {
            if (pre == 's' || pre == 'c')
            {
               if (vlen == 16)
                  VECi = VTAVXZ;
               else if (vlen == 8)
                  VECi = VTAVX;
               else if (vlen == 4 && VECi == VTAVX)
                  VECi = VTSSE;
               else if (vlen == 4)
                  VECi = VTGV;
               else
                  VECi = VTSC;
            }
            else
            {
               if (vlen == 8)
                  VECi = VTAVXZ;
               else if (vlen == 4)
                  VECi = VTAVX;
               else if (vlen == 2 && VECi == VTAVX)
                  VECi = VTSSE;
               else if (vlen == 2)
                  VECi = VTGV;
               else
                  VECi = VTSC;
            }
         }
         mp->genstr = GenString(pre, mp->lat, mp->kbB, mu, mp->nu, mp->ku,
                                mp->kmaj, mp->rout);
         VECi = vl0;
      }
   }
}

ATL_mmnode_t *GetNewGenNode(char pre, int nb, int lat, int mu, int nu, int ku,
                            int kmaj)
{
   ATL_mmnode_t *np;
   char *ln;
   if (kmaj < 2)
      kmaj = 0;
   np = GetMMNode();
   np->rout = malloc(sizeof(char)*27);
   assert(np->rout);
   sprintf(np->rout, "ATL_%cgamm%d_%dx%d_nb%d.c", pre, kmaj, mu, nu, nb);
   np->vlen = VLEN[VECi];
   np->nbB = np->mbB = np->kbB = nb;
   if (kmaj)
      np->mu = mu;
   else
      np->mu = mu*VLEN[VECi];
   np->nu = nu;
   np->ku = ku;
   np->muladd = VECi;  /* stores what vector ISA to use */
   np->lat = lat;
   np->kmaj = kmaj;
   np->genstr = GenString(pre, lat, nb, mu, nu, ku, kmaj, np->rout);
   ApplyMoves2Flag(np, IMVS);
   return(np);
}

ATL_mmnode_t *FindDefMUNU(int verb, char pre, int nreg, int lat, int nb, int ku,
                          int *MU, int *NU)
{
   ATL_mmnode_t *mmp;
   double mf, mfB=0.0;
   int n, i, j, kb, muB=1, nuB=1;

   mmp = ReadMMFileWithPath(pre, "res", "gAMMMUNU.sum");
   if (mmp)
   {
      FillInGenStrings(pre, mmp);
      nb = mmp->kbB;
      if (mmp->mflop[0] < 0.0)
         mmp->mflop[0] = TimeMMKernel(verb, 1, mmp, pre, nb, nb, nb,
                                      nb, nb, nb, 1, 0, -1);
      printf("READ IN BEST GENNED MU=%d, NU=%d, MFLOP=%.2f\n\n",
             mmp->mu, mmp->nu, mmp->mflop[0]);
/*
 *    See if there is a mismatch between vector settings
 */
      if (mmp->vlen != VLEN[VECi])
      {
         printf("\n\n!!! WARNING: TURNING OFF VECTORIZATION DUE TO MISMATCHED VLEN IN 'res/%cAMMMUNU.sum!!!!\n\n", pre);
         VECi = VTSC;
      }
      *MU = mmp->mu / VLEN[VECi];
      assert(*MU);
      *NU = mmp->nu;
      return (mmp);
   }
   mmp = GetMMNode();
   mmp->ku = ku;
   mmp->muladd = (lat != 0);
   mmp->lat = lat;
   mmp->rout = DupString("ATL_Xamm_munu.c");
   mmp->rout[4] = pre;
   mmp->mbB = mmp->nbB = mmp->kbB = nb;
   mmp->vlen = VLEN[VECi];
/*
 * Try all near-square register blocking cases
 */
   printf("Finding best MUxNU case for nb=%d\n", nb);
   for (n=4; n < nreg; n++)
   {
      int mbu, nbu, mu, nu;
      for (j=1; j*j < n; j++);
      i = n / j;
      if (nb%i || nb%j)
         continue;
      mu = mmp->mu = i * VLEN[VECi];
      nu = mmp->nu = j;
      if (mmp->genstr)
        free(mmp->genstr);
      mbu = (nb >= mu) ? (nb/mu)*mu : mu;
      nbu = (nb >= nu) ? (nb/nu)*nu : nu;
      mmp->genstr = GenString(pre, lat, nb, i, j, ku, mmp->kmaj, mmp->rout);
      mf = TimeMMKernel(verb, 1, mmp, pre, mbu, nbu, nb, nb, nb, nb, 1, 0, -1);
      printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
      if (mf > mfB)
      {
         muB = i;
         nuB = j;
         mfB = mf;
      }
   }
/*
 * For x86, try 1-D cases since older machines are 2-operand assemblies
 */
   #if defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)
   if (VECi != VTAVX && VECi != VTAVXZ)  /* AVX is 3-operand */
   {
      printf("BEST NEAR-SQUARE CASE IS MU=%d, NU=%d, MFLOP=%.2f\n\n",
             muB, nuB, mfB);
      printf("Finding best 1-D outer loop unrolling for nb=%d\n", nb);
      for (n=2; n < nreg; n++)
      {
         int mbu, nbu, mu, nu;
         i = 1; j = n;
         if (nb % n)
            continue;
         mu = mmp->mu = i*VLEN[VECi];
         nu = mmp->nu = j;
         if (mmp->genstr)
           free(mmp->genstr);
         mmp->genstr = GenString(pre, lat, nb, i, j, ku, mmp->kmaj, mmp->rout);
         mbu = (nb >= mu) ? (nb/mu)*mu : mu;
         nbu = (nb >= nu) ? (nb/nu)*nu : nu;
         mf = TimeMMKernel(verb, 1, mmp, pre, mbu, nbu, nb, nb, nb, nb,
                           1, 0, -1);
         printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
         if (mf > mfB)
         {
            muB = i;
            nuB = j;
            mfB = mf;
         }
         i = n; j = 1;
         mu = mmp->mu = i * VLEN[VECi];
         nu = mmp->nu = j;
         mbu = (nb >= mu) ? (nb/mu)*mu : mu;
         nbu = (nb >= nu) ? (nb/nu)*nu : nu;
         if (mmp->genstr)
           free(mmp->genstr);
         mmp->genstr = GenString(pre, lat, nb, i, j, ku, mmp->kmaj, mmp->rout);
         mf = TimeMMKernel(verb, 1, mmp, pre, mbu, nbu, nb, nb, nb, nb,
                           1, 0, -1);
         printf("   MU=%2d, NU=%2d, MFLOP=%.2f\n", i, j, mf);
         if (mf > mfB)
         {
            muB = i;
            nuB = j;
            mfB = mf;
         }
      }
   }
   #endif
   mmp->mu = muB * VLEN[VECi];
   mmp->nu = nuB;
   mmp->mbB = (nb >= muB) ? (nb/muB)*muB : muB;
   mmp->nbB = (nb >= nuB) ? (nb/nuB)*nuB : nuB;
   mmp->kbB = nb;
   if (mmp->genstr)
     free(mmp->genstr);
   mmp->genstr = GenString(pre, lat, nb, muB, nuB, ku, mmp->kmaj, mmp->rout);
   WriteMMFileWithPath(pre, "res", "gAMMMUNU.sum", mmp);
   printf("BEST CASE IS MU=%d, NU=%d, MFLOP=%.2f (%.2f)\n\n",
          muB, nuB, mf, mfB);
   *MU = muB;
   *NU = nuB;
   return(mmp);
}

void GetMUNUbyNB(int nb, int nreg, int *MU, int *NU)
{
   int mu=(*MU), nu=(*NU), vmu=mu*VLEN[VECi];

   assert(mu && nu && !(nb%VLEN[VECi]));
   if (!(nb%vmu) && !(nb%nu))
      return;
   if (nu == 1) /* handle MUx1 by decreasing by VLEN */
   {
      int u = vmu, vlen = VLEN[VECi];
      while (u+u+1 <= nreg && nb%u)
         u += vlen;
      if (u+u+1 > nreg)
         u -= vlen;
      while (nb%u)
         u -= vlen;
      assert(u);
      *MU = u / vlen;
      return;
   }
   if (mu == 1 || nu == 1) /* handle 1-D cases by just inc/dec U */
   {
      int u = (mu == 1) ? nu : mu;
      while (u+u+1 <= nreg && nb%u)
         u++;
      if (u+u+1 > nreg)
         u--;
      while (nb%u)
         u--;
      if (mu == 1)
         *NU = u;
      else
         *MU = u;
      return;
   }
   if (nb%vmu)  /* mu can't handle NB */
   {
      int i;
/*
 *    try increasing mu until we run out of registers
 */
      for (i=mu+1; i*nu+i+nu <= nreg; i++)
         if (!(nb%(i*VLEN[VECi])))
            break;
/*
 *    Try decreasing mu until it divides
 */
      if (nb%(i*VLEN[VECi]) || i*nu+i+nu > nreg)
      {
         for (mu--; mu; mu--)
            if (!(nb%(mu*VLEN[VECi])))
               break;
      }
      else
         mu = i;
   }
   if (nb%nu)  /* nu can't handle NB */
   {
      int i;
/*
 *    try increasing nu until we run out of registers
 */
      for (i=nu+1; i*mu+i+mu <= nreg; i++)
         if (!(nb%i))
            break;
/*
 *    Try decreasing nu until it divides
 */
      if (nb%i || i*mu+i+mu > nreg)
      {
         for (nu--; nu; nu--)
            if (!(nb%nu))
               break;
      }
      else
         nu = i;
   }
   *MU = mu;
   *NU = nu;
}

int FindNBInArray(int nb, int *nbs)
/*
 * RETURNS: location+1, or 0 if not found
 */
{
   int i, n = (nbs[0] > 0) ? nbs[0] : -nbs[0];
   for (i=1; i <= n; i++)
       if (nbs[i] == nb)
          return(i);
   return(0);
}
ATL_mmnode_t *CreateGenCasesFromNBs
(
   ATL_mmnode_t *mmb,   /* best user-contributed kernels */
   char pre,            /* precision: s/d */
   int *nbs,            /* list of desired NBs */
   int nreg,            /* upper bound on register use */
   int MU, int NU,      /* default M/N unrolling */
   int KU               /* -1 for fully unrolled, else unrolling factor */
)
/*
 * Generate a list of generated kernels, with the union of nb's in nbs
 * and mmb, and return the generated nodes for timing.
 * HERE HERE HERE: this code is crap, needs to merge both lists, not user
 * list twice.
 */
{
   ATL_mmnode_t *mp, *umb=NULL, *ap;
   int i, n = (nbs[0] > 0) ? nbs[0] : -nbs[0], ne=0, *enbs;

/*
 * Create new queue with an entry for all NBs; both lists (mmb & nbs) are
 * sorted in increasing size
 */
   if (!n && !mmb)
      return(NULL);
   n++;
   ap = mmb;  /* add ptr */
   i = 1;     /* ptr to normal block under consideration */
   do
   {
      int nb, mu=MU, nu=NU, ku;
      ATL_mmnode_t *p=NULL;
      if (ap && i < n)  /* must choose amongst blocks */
      {
         nb = ap->kbB;
         nb = Mmin(nb, nbs[i]);
         if (nb == ap->kbB)
            ap = ap->next;
         if (nb == nbs[i])
            i++;
      }
      else if (ap)
      {
         nb = ap->kbB;
         ap = ap->next;
      }
      else
         nb = nbs[i++];
      ku = (KU == -1) ? nb : KU;

/*
 *    If NB is not a multiple of VLEN, drop down to shorter ops
 */
      if (nb%VLEN[VECi])
      {
         int vl=VECi;
/*
 *       For AVX, see if dropping to SSE will fix problem
 */
         if (VECi == VTAVX && !(nb%VLEN[VTSSE]))  /* AVX can drop to SSE */
         {
            VECi = VTSSE;
            GetMUNUbyNB(nb, nreg, &mu, &nu);
            p = GetNewGenNode(pre, nb, 0, mu, nu, ku, 0);
         }
         if (!p)
         {
            VECi = VTSC;
            GetMUNUbyNB(nb, nreg, &mu, &nu);
            p = GetNewGenNode(pre, nb, 0, mu, nu, ku, 0);
         }
         VECi = vl;
      }
      else
      {
         GetMUNUbyNB(nb, nreg, &mu, &nu);
         p = GetNewGenNode(pre, nb, 0, mu, nu, ku, 0);
      }
      if (umb)
      {
         mp->next = p;
         mp = p;
      }
      else
         umb = mp = p;
   }
   while (i < n || ap);
   return(umb);
}
void SetGenVec(int verb, char pre)
/*
 * This routine uses a simple timing to be sure if vectorization helps or not
 */
{
   ATL_mmnode_t *mp;
/*
 * If vector operations are being used, make sure they work; compiler and
 * flag changes can mess them up, and in this case we'll fall back to
 * scalar generation.  Try to see if we can successfully test simplist
 * possible vector kernel, and fall back to scalar kernel if we can't
 */
   if (VLEN[VECi] < 2)
      return;
   mp = GetNewGenNode(pre, 32, 0, 1, 1, 1, 0);
   if (MMKernelFailsTest(pre, 32, 32, 32, 1, mp))
   {
      printf("ERROR: VEC='%s' FAILED, genstr='%s'!\n",VECs[VECi],mp->genstr);
      KillMMNode(mp);
/*
 *    For AVX, try falling back to SSE
 */
      if (VECi == VTAVX)
      {
         VECi = VTSSE;
         KillMMNode(mp);
         mp = GetNewGenNode(pre, 32, 0, 1, 1, 1, 0);
         if (MMKernelFailsTest(pre, 32, 32, 32, 1, mp))
            VECi = VTSC;
      }
      else
         VECi = VTSC;
   }
   KillMMNode(mp);
/*
 * For AVX, switch to SSE if AVX doesn't offer a performance advantage
 * (as on AMD Dozer), since SSE smaller code size and requires less cleanup
 */
   if (VECi == VTAVX)
   {
      double mfA, mfS, mf;
      char *sp;
      int vl;
      mp = GetNewGenNode(pre, 128, 0, 1, 4, 1, 0);
      mfA = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 128, 128, 128,
                         1, 0, -1);
      KillMMNode(mp);
      mp = GetNewGenNode(pre, 128, 0, 2, 2, 1, 0);
      mf = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 128, 128, 128,
                        1, 0, -1);
      KillMMNode(mp);
      if (mf > mfA)
         mfA = mf;
      vl = VECi;
      VECi = VTSSE;
      mp = GetNewGenNode(pre, 128, 0, 1, 4, 1, 0);
      mfS = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 128, 128, 128,
                         1, 0, -1);
      KillMMNode(mp);
      mp = GetNewGenNode(pre, 128, 0, 2, 2, 1, 0);
      mf = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 128, 128, 128,
                        1, 0, -1);
      if (mf > mfA)
         mfA = mf;
      KillMMNode(mp);
      if (mfA < 1.03*mfS)
         printf("USING SSE INSTEAD OF AVX, AVX=%.2f, SSE=%.2f\n", mfA, mfS);
      else
      {
         printf("AVX GOOD, AVX=%.2f, SSE=%.2f\n", mfA, mfS);
         VECi = vl;
      }
   }
/*
 * For any system, don't use vector instructions if they aren't faster than
 * scalar.
 */
   if (VLEN[VECi] > 1)
   {
      double mfV, mfS;
      char *sp;
      int vl;
      mp = GetNewGenNode(pre, 128, 0, 1, 4, 1, 0);
      mfV = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 128, 128, 128,
                         1, 0, -1);
      KillMMNode(mp);
      vl = VECi;
      VECi = VTSC;
      mp = GetNewGenNode(pre, 128, 0, 1, 4, 1, 0);
      mfS = TimeMMKernel(verb, 1, mp, pre, 128, 128, 128, 128, 128, 128,
                         1, 0, -1);
      KillMMNode(mp);
      if (mfV < 1.05*mfS)
         printf("USING SCALAR INSTEAD OF VECTOR, VEC=%.2f, SCALAR=%.2f\n",
                mfV, mfS);
      else
      {
         printf("VEC GOOD, VEC=%.2f, SCALAR=%.2f\n", mfV, mfS);
         VECi = vl;
      }
   }
   printf("GENERATING WITH VEC='%s', VLEN=%d\n\n", VECs[VECi], VLEN[VECi]);
}

ATL_mmnode_t *FindBestGenCases(int verb, char pre, int nreg,
                               int *nbs, ATL_mmnode_t *ummb)
{
   ATL_mmnode_t *mp, *gmmU, *gmmb;
   int MU, NU;

   gmmb = ReadMMFileWithPath(pre, "res", "gAMMRES.sum");
   if (gmmb)
   {
      printf("Reading in generated cases for all NBs:\n");
      FillInGenStrings(pre, gmmb);
      for (mp=gmmb; mp; mp = mp->next)
      {
         const int mu = (mp->vlen) ? mp->mu / mp->vlen : mp->mu;
         int kb = mp->kbB, mb = Mmax(mp->mu,kb), nb=Mmax(mp->nu,kb);
         if (mp->mflop[0] < 0.0)
            mp->mflop[0] = TimeMMKernel(verb, 1, mp, pre, mb, nb, kb,
                                         kb, kb, mb, 1, 0, -1);
         printf("  NB=%d, MU=%d, NU=%d, vlen=%d, MFLOP=%.2f\n",
                nb, mu, mp->nu, mp->vlen, mp->mflop[0]);
      }
      WriteMMFileWithPath(pre, "res", "gAMMRES.sum", gmmb);
      printf("Done.\n\n");
      return(gmmb);
   }
   SetGenVec(verb, pre);
/*
 * Find the best mu/nu for NB=120; we don't care if we overflow cache for
 * this timing, and 120 = LCM(2,3,4,5,6,8,12).  Use ku=1 so that large
 * problems don't have large K-driven advantage.
 */
   KillMMNode(FindDefMUNU(verb, pre, nreg, 0, 120, 1, &MU, &NU));
   gmmb = CreateGenCasesFromNBs(ummb, pre, nbs, nreg, MU, NU, 1);
   if (verb > 1)
   {
      printf("\n");
      PrintMMNodes(stdout, gmmb);
      printf("\n");
   }
   printf("Finding generated cases for all NBs:\n");
   for (mp=gmmb; mp; mp = mp->next)
   {
      int mu;
      mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mp->mbB, mp->nbB, mp->kbB,
                                  mp->kbB, mp->kbB, mp->mbB, 1, 0, -1);
      mu = (mp->vlen) ? mp->mu / mp->vlen : mp->mu;
      printf("   NB=%d, mu=%d, nu=%d, vlen=%d, MFLOPS=%.2f\n",
             mp->kbB, mu, mp->nu, mp->vlen, mp->mflop[0]);
   }
   printf("Done.\n\n");
   WriteMMFileWithPath(pre, "res", "gAMMRES.sum", gmmb);
   return(gmmb);
}


ATL_mmnode_t *GetWorkingUserCases(int verb, char pre)
{
   ATL_mmnode_t *mmb, *mmp;
   mmb = ReadMMFileWithPath(pre, "res", "WORKING.sum");
   if (mmb)
      return(mmb);
   mmb = ReadMMFileWithPath(pre, "AMMCASES", "amcases.idx");
   if (!mmb)
      return(mmb);
/*
 * Eliminate those kernels that can't work for any block size
 */
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      if (FLAG_IS_SET(mmp->flag, MMF_KUISKB))
         mmp->mbB = mmp->nbB = mmp->kbB = mmp->ku;
      else
      {
         int m = Mylcm(mmp->mu, mmp->nu);
         m = ((60+m-1)/m)*m;
         mmp->mbB = mmp->nbB = mmp->kbB = m;
         if (mmp->kbmin)
            mmp->kbB = Mmax(mmp->kbB, mmp->kbmin);
         if (mmp->kbmax)
            mmp->kbB = Mmin(mmp->kbB, mmp->kbmax);
      }
   }
   mmb = DelBadMMKernels(pre, verb, mmb);
   WriteMMFileWithPath(pre, "res", "WORKING.sum", mmb);
   return(mmb);
}

ATL_mmnode_t *FindBestUserCases(int verb, char pre, int *nbs, ATL_mmnode_t *mmb)
/*
 * NOTE: frees mmb after search!!
 * RETURNS: list of the best user case for each supplied NB; if no user case
 *          works, special "generated" node is returned for later filling out.
 */
{
   ATL_mmnode_t *mmp, *mp;
   mp = ReadMMFileWithPath(pre, "res", "uAMMRES.sum");
/*
 * If final output file exists, then we need to rerun timings at worst
 */
   if (mp)
   {
      KillAllMMNodes(mmb);
      for (mmp=mp; mmp; mmp = mmp->next)
      {
         if (mmp->ID > 0 && mmp->mflop[0] < 0.0)
            mmp->mflop[0] = TimeMMKernel(verb, 0, mmp, pre, mmp->mbB, mmp->nbB,
                                         mmp->kbB, mmp->kbB, mmp->kbB, mmp->mbB,
                                         1, 0, -1);
         if (mmp->ID > 0)
            printf("USER KERNEL AT NB=%d gets MFLOP=%.2f\n",
                   mmp->kbB, mmp->mflop[0]);
         else
            printf("NO USER KERNEL FOR NB=%d\n", mmp->kbB);
      }
      printf("\n");
      return(mp);
   }
   mmp = FindBestForEachNB(verb, pre, mmb, nbs);
   KillAllMMNodes(mmb);
   WriteMMFileWithPath(pre, "res", "uAMMRES.sum", mmp);
   return(mmp);
}

ATL_mmnode_t *MergeCases
(
   int imf,
   ATL_mmnode_t *bs0, /* queue of cases */
   ATL_mmnode_t *bs1  /* queue of cases */
)
/*
 * Merges two queues of matmul kern cases.  Cases are not winnowed, but
 * duplicates are not allowed, so if two entries have the same kbB, then
 * we take the one with best mflop[imf].  If imf < 0, then we do indeed
 * allow duplicates of kbB.
 * NOTE: does not change bs0 or bs1.
 * ASSUMES: both bs0 & bs1 are in kb-increasing order.
 * RETURNS: base ptr to merged queue
 */
{
   ATL_mmnode_t *mb=NULL, *mp;
   while (bs0 || bs1)
   {
      ATL_mmnode_t *p;
      if (bs0 && bs1)
      {
         if (bs0->kbB < bs1->kbB)
         {
            p = CloneMMNode(bs0);
            bs0 = bs0->next;
         }
         else if (bs0->kbB > bs1->kbB)
         {
            p = CloneMMNode(bs1);
            bs1 = bs1->next;
         }
         else /* they are equal, must take best performer, or both */
         {
/*
 *          If we are taking both, special case can't use general completion
 */
            if (imf < 0)
            {
               p = CloneMMNode(bs0);
               bs0 = bs0->next;
               p->next = CloneMMNode(bs1);
               bs1 = bs1->next;
               if (mb)
                  mp->next = p;
               else
                  mb = p;
               mp = p->next;
               continue;
            }
/*
 *          Taking only the best performer, but moving both base ptrs
 */
            else
            {
/*
 *             If they are equal, take the KRUN=1 case if it exists, else
 *             take the most flexible one or one requiring the least cleanup
 */
               if (bs0->mflop[imf] == bs1->mflop[imf])
               {
                  if (FLAG_IS_SET(bs0->flag, MMF_KRUNTIME))
                     p = bs0;
                  else if (FLAG_IS_SET(bs1->flag, MMF_KRUNTIME))
                     p = bs1;
                  else if (bs0->ku < bs1->ku)
                     p = bs0;
                  else if (bs1->ku < bs0->ku)
                     p = bs1;
                  else
                  {
                     const int u0=Mmax(bs0->mu, bs0->nu),
                               u1=Mmax(bs1->mu, bs1->nu);
                     p = (u0 <= u1) ? bs0 : bs1;
                  }
               }
               else
                  p = (bs0->mflop[imf] > bs1->mflop[imf]) ? bs0 : bs1;
               p = CloneMMNode(p);
               bs0 = bs0->next;
               bs1 = bs1->next;
            }
         }
      }
      else if (bs0)
      {
         p = CloneMMNode(bs0);
         bs0 = bs0->next;
      }
      else /* if (bs1) */
      {
         p = CloneMMNode(bs1);
         bs1 = bs1->next;
      }
      if (mb)
      {
         mp->next = p;
         mp = p;
      }
      else
        mp = mb = p;
   }
   return(mb);
}

#define HUGE_NB 180
ATL_mmnode_t *WinnowHugeNB
(
   int imf,
   ATL_mmnode_t *mb  /* queue of cases */
)
/*
 * Removes any NB >= HUGE_NB that aren't at least 2% faster than smaller cases
 */
{
   ATL_mmnode_t *mp, *p, *prev=mb;
   double mfB;

   if (!mb || !mb->next)
      return(mb);
   mp = mb->next;
/*
 * Find best-performing kernel below HUGE_NB
 */
   mfB = mp->mflop[imf];
   for (mp=mb->next; mp; mp = mp->next)
   {
      if (mp->mbB < HUGE_NB && mp->nbB < HUGE_NB && mp->kbB < HUGE_NB)
         mfB = Mmax(mfB, mp->mflop[imf]);
      else
         break;
      prev = mp;
   }
/*
 * If no kernels above threshold, return original queue
 */
   if (!mp)
      return(mb);
/*
 * mp points to first NB above threshold, but there is no point in deleting
 * small NB if we leave large NB, so delete only from the end of queue
 */
  do
  {
     for (p=mp; p->next; p = p->next);
     if (p->mflop[imf] <= 1.02*mfB)
        mp = RemoveMMNodeFromQ(mp, p);
     else  /* stop removing stuff */
        break;
  }
  while (mp);
  prev->next = mp;
  return(mb);
}

ATL_mmnode_t *WinnowCases
(
   int imf,
   ATL_mmnode_t *mb  /* queue of cases */
)
/*
 * Removes any case that runs slower than a smaller case
 * RETURNS: mb with queue bad kernels deleted
 * NOTE: mb can never change, since by def nothing smaller than 1st case
 */
{
   ATL_mmnode_t *prev = mb, *mp;

   if (!mb)
      return(NULL);
   mp = mb->next;
   while (mp)
   {
      if (mp->mflop[imf] <= prev->mflop[imf])  /* kill slow KB */
         mp = prev->next = KillMMNode(mp);
      else
      {
         prev = mp;
         mp = mp->next;
      }
   }
   return(mb);
}

ATL_mmnode_t *MergeAndWinnowCases
(
   int verb,
   char pre,
   ATL_mmnode_t *umb, /* queue of user cases */
   ATL_mmnode_t *gmb  /* genned cases, always include NBs of umb */
)
/*
 * Merges user and gmp cases, while getting rid of cases that get worse
 * performance than their smaller blocks; FREES umb and gmb
 * RETURNS: new merged and winnowed queue
 */
{
   ATL_mmnode_t *mmb=NULL, *mmp, *gmp, *ump=umb;
   for (gmp=gmb; gmp; gmp = gmp->next)
   {
      ATL_mmnode_t *p;
      if (ump)
      {
         if (ump->kbB == gmp->kbB)
         {
            if (gmp->mflop[0] >= ump->mflop[0])
               p = CloneMMNode(gmp);
            else
               p = CloneMMNode(ump);
            ump = ump->next;
         }
         else
            p = CloneMMNode(gmp);
      }
      else
         p = CloneMMNode(gmp);
      p->next = NULL;
      if (mmb)
      {
/*
 *       If larger NB isn't faster than smaller one, kill it for nb >= 16
 */
         if (p->kbB >= 16 && mmp->mflop[0] > p->mflop[0])
            KillMMNode(p);
         else
         {
            mmp->next = p;
            mmp = p;
         }
      }
      else
         mmp = mmb = p;
   }
   KillAllMMNodes(umb);
   KillAllMMNodes(gmb);
   mmb = WinnowCases(0, mmb);
   return(mmb);
}



int FailKCleanTests(char pre, int nb, ATL_mmnode_t *kp)
/*
 *  This routine tests if a kernel is suitable for use in K-cleanup by
 *  doing testing with ku=1, kb=0, and tries all K values between 1 and nb
 *  RETURNS: 0 if kernel passes all tests, else non-zero
 */
{
   int i, beg, end, inc;

   if (!FLAG_IS_SET(kp->flag, MMF_KRUNTIME) ||
       (kp->ku != 1 && kp->ku != kp->kmaj))
      return(-1);
   printf("TESTING ID=%d, rout='%s', nb=%d, mu=%d, nu=%d for K-cleanup:\n",
          kp->ID, kp->rout, nb, kp->mu, kp->nu);

   if (kp->kmaj > 1)
   {
      inc = beg = kp->kmaj;
      end = ((nb+inc-1)/inc)*inc;
   }
   else
   {
      beg = inc = 1;
      end = nb;
   }
   for (i=beg; i <= end; i += inc)
   {
      int ierr;
      ierr = MMKernelFailsTest(pre, nb, nb, i, 0, kp);
      if (ierr)
      {
         printf("  K=%d: FAILED!\n", i);
         return(ierr);
      }
      else
         printf("  K=%d: PASSED!\n", i);
   }
   printf("PASSED ALL K-tests!\n\n");
   return(0);
}
ATL_mmnode_t *GetUniqueKClean(int verb, char pre, ATL_mmnode_t *mmb)
/*
 * OUTPUT: <pre>UMMKCLEAN.sum: all unique kerns to be compiled
 */
{
   ATL_mmnode_t *mp, *gmmb, *ummb, *ub, *np, **dlmm;
   int nn=0, nd=0, n=0;  /* #needed & done, total, copy of done */
   int *dl, *nl;         /* done and needed lists */
   int i;
/*
 * Find out how many total kernels, and how many already have their own
 * cleanup (nd, number done).  This nd may be bigger than it should, because
 * we can't guarantee they are unique
 */
   for (mp=mmb; mp; mp = mp->next, n++)
      if ((mp->ku == 1 || (mp->kmaj == mp->ku)) &&
          FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
        nd++;
   dl = malloc(8*n*sizeof(int));
   assert(dl);
   if (nd)
   {
      dlmm = malloc(nd*sizeof(ATL_mmnode_t*));
      assert(dlmm);
   }
   else
      dlmm = NULL;

   nl = dl + (n<<2);
   nd = 0;
/*
 * First, go back through kernels, and add kernels that can serve as K-cleaners
 * to the done list
 */
   for (mp=mmb; mp; mp = mp->next, n++)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME) &&
          (mp->ku == 1 || mp->kmaj == mp->ku))
      {
         const int nd4 = (nd<<2), mu=mp->mu, nu=mp->nu, kmaj=mp->kmaj;
/*
 *       See if trip is already in done list if so, no new entry, just update kb
 *       and cleanup kernel entry
 */
         for (i=0; i < nd4; i += 4)
            if (mu == dl[i] && nu == dl[i+1] && kmaj == dl[i+2])
               break;
         if (i < nd4)
         {                      /* (larger NB always later) */
            dl[i+3] = mp->kbB;  /* take largest kbB that matches mu/nu */
            dlmm[i>>2] = mp;
            continue;
         }
         else
         {
            dl[nd4] = mu;
            dl[nd4+1] = nu;
            dl[nd4+2] = kmaj;
            dl[nd4+3] = mp->kbB;
            dlmm[nd++] = mp;
         }
      }
   }
/*
 * Delete any kernels from dl that fail to actually work for K cleaning
 */
   for (i=0; i < nd; i++)
   {
      if (FailKCleanTests(pre, dlmm[i]->kbB, dlmm[i]))
      {
         const int i4=(i<<2), nc=nd-i-1;
         if (nc > 0)
         {
            memcpy(dl+i4, dl+i4+4, (nc<<2)*sizeof(int));
            memcpy(dlmm[i], dlmm[i+1], nc*sizeof(ATL_mmnode_t*));
         }
         nd--;
      }
   }

/*
 * Find all unique (mu,nu,kmaj) combos that still need to to be cleaned;
 * there will be nn (# needed) of these, and we'll save (mu,nu,MAXNB) in
 * needed list (nl).
 * We use MAXNB for testing (large NB tests mosts cases of K).
 * Combos that are handled by the done list (dl) aren't added to needed list.
 */
   for (mp=mmb; mp; mp = mp->next, n++)
   {
      int mu=mp->mu, nu=mp->nu, kmaj=mp->kmaj, nn4=(nn<<2), nd4=(nd<<2);
/*
 *    See if pair is already in done list or needed list, if so, no change
 */
      for (i=0; i < nd4; i += 4)
         if (mu == dl[i] && nu == dl[i+1] && kmaj == dl[i+2])
            break;
      if (i < nd4)    /* if it was found in the done list */
         continue;    /* this combo is already handled */
/*
 *    If we reach here, combo is not handled, must add to needed list
 */
      for (i=0; i < nn4; i += 4)
         if (mu == nl[i] && nu == nl[i+1] && kmaj == nl[i+2])
            break;
      if (i < nn4)            /* If already in needed list */
      {
         nl[i+3] = mp->kbB;   /* just update kb so we get largest for testing */
         continue;
      }
/*
 *    If we haven't seen this pair before, add to needed list
 */
      else
      {
         nl[nn4] = mu;
         nl[nn4+1] = nu;
         nl[nn4+2] = kmaj;
         nl[nn4+3] = mp->kbB;
         nn++;
      }
   }
/*
 * Now, create a queue of generated kernels for each needed pair, and time
 * it's maxNB performance.
 */
   gmmb = NULL;
   printf("Timing Generated K-cleanup:\n");
   for (i=0; i < nn; i++)
   {
      ATL_mmnode_t *p;
      const int i4 = (i<<2), mu=nl[i4], nu=nl[i4+1], kmaj=nl[i4+2];
      int nb = Mmax(nl[i4+3],nu), mb = (nb > mu) ? (nb/mu)*mu : mu;
      const int kb = (nb > 8) ? (nb>>2) : nb;
      int vl=VECi, vmu, KK;
      double mf;
/*
 *    HERE HERE: Improve KMAJ when generator is extended!
 */
      if (kmaj > 1)
      {
        vmu = mu;
        VECi = VTSC;
      }
      else
      {
         if (mu % VLEN[VECi])
         {
            if (VECi == VTAVX && !(mu%VLEN[VTSSE]))
               VECi = VTSSE;
            else
               VECi = VTSC;
         }
         vmu = mu / VLEN[VECi];
      }
/* HERE HERE */
      p = GetNewGenNode(pre, nb, 0, vmu, nu, 1, kmaj);
      p->mbB = mb;
      p->nbB = nb;
      p->kbB = kb;
      p->flag |= (1<<MMF_KRUNTIME);
      VECi = vl;
      #if 0  /* by default don't waste time testing generated code */
         assert(!FailKCleanTests(pre, nb, p));
      #endif
      KK = (kmaj < 2) ? kb : ((kb+kmaj-1)/kmaj)*kmaj;
      p->mflop[0] = TimeMMKernel(verb, 0, p, pre, mb, nb, kb, 0, 0, nb,
                                 0, 0, -1);
      if (KK != kb)
         p->mflop[0] *= (double)kb / (double)KK;
      printf("   nb=%d, kb=%d,  mu=%d, nu=%d, MFLOP=%.2f\n",
             nb, kb, mu, nu, p->mflop[0]);
      if (gmmb)
      {
         mp->next = p;
         mp = p;
      }
      else
         gmmb = mp = p;
   }
   printf("Done.\n");
/*
 * Now, add the done-list items to generated list
 */
   for (i=0; i < nd; i++)
   {
      const int i4=4*i, mu=dl[i4], nu=dl[i4+1], kmaj=dl[i4+2], nb=dl[i4+3];
      ATL_mmnode_t *prev=NULL;
/*
 *    Get a copy of done-list kern that can be added to genlist
 */
      np = CloneMMNode(dlmm[i]);
      np->next = gmmb;
      gmmb = np;
   }
   if (dlmm)
      free(dlmm);
/*
 * Now, search index file for suitable user-submitted kernels to compete
 * with existing solutions
 */
   ub = ReadMMFileWithPath(pre, "AMMCASES", "amcases.idx");
   ummb = NULL;  /* no suitable user cases to begin */
/*
 * Look through user-list for any routine with ku=1 and K-Runtime
 */
   for (mp=ub; mp; mp = mp->next)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME) &&
          (mp->ku == 1 || mp->kmaj == mp->ku))
      {
/*
 *       It matched our gross criteria, see if it is a required mu/nu
 */
         for (i=0; i < nn; i++)
         {
            const int i4=(i<<2), mu=nl[i4], nu=nl[i4+1], kmaj=nl[i4+2],
                      nb=nl[i4+3];
            if (mp->mu == mu && mp->nu == nu && mp->kmaj == kmaj)
            {
               if (!FailKCleanTests(pre, nb, mp))
               {
                  ATL_mmnode_t *p;
                  p = CloneMMNode(mp);
                  p->next = NULL;
                  p->nbB = ((nb+nu-1)/nu)*nu;
                  p->mbB = ((nb+mu-1)/mu)*mu;
                  p->kbB = nb;
                  if (ummb)
                  {
                     np->next = p;
                     np = p;
                  }
                  else
                     np = ummb = p;
                  break;
               }
            }
         }
      }
   }
   KillAllMMNodes(ub);
/*
 * If we have both user and genned code, must compare timing to select best
 */
   if (ummb)
   {
/*
 *    Now, loop over user cases and time them for comparison with genned
 */
      printf("Timing User K-cleanup:\n");
      for (mp=ummb; mp; mp = mp->next)
      {
         const int nb = mp->nbB, kb = (nb > 8) ? (nb>>2) : nb;
         const int KK = (mp->kmaj < 2) ? kb:((kb+mp->kmaj-1)/mp->kmaj)*mp->kmaj;
         mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mp->mbB, nb, KK, 0, 0,
                                     KK, 0, 0, -1);
         if (KK != kb)
            mp->mflop[0] *= (double)kb / (double)KK;
         printf("   ID=%d, nb=%d, kb=%d, mu=%d, nu=%d, MFLOP=%.2f\n",
                mp->ID, nb, kb, mp->mu, mp->nu, mp->mflop[0]);
      }
      printf("Done timing, merging lists:\n");
/*
 *    Merge generated (gmmb) and user (ummb) kerns by selecting best performing.
 *    gmmb is a superset of ummb, so what we will do is look through gmmb
 *    for matching (mu,nu,dup), time them, and if ummb is faster, replace
 *    that entry in gmmb with ummb.
 */
      while (ummb)
      {
         ATL_mmnode_t *prev=NULL;
         int mu=ummb->mu, nu=ummb->nu, kmaj = ummb->kmaj;
         for (mp=gmmb; mp && (mp->mu != mu || mp->nu != nu || mp->kmaj != kmaj);
              mp = mp->next)
            prev = mp;
         assert(mp);  /* logic error if we haven't found it */
/*
 *       If user case gets better performance, replace genned case in queue
 */
         if (ummb->mflop[0] > gmmb->mflop[0])
         {
            printf("   Replacing genned case (%.2f) with user ID %d (%.2f)\n",
                   gmmb->mflop[0], ummb->ID, ummb->mflop[0]);
            if (prev)
            {
               prev->next = ummb;
               ummb = ummb->next;
               prev->next->next = KillMMNode(mp);
            }
            else /* replace gmmb, mp pts at gmmb */
            {
               ATL_mmnode_t *up=ummb;
               ummb = ummb->next;
               up->next = KillMMNode(gmmb);
               gmmb = up;
            }
         }
         else /* user case loser, just delete it */
         {
            printf("   Preferring genned case (%.2f) over user ID %d (%.2f)\n",
                   gmmb->mflop[0], ummb->ID, ummb->mflop[0]);
            ummb = KillMMNode(ummb);
         }
      }
      printf("DONE.\n\n");
   }
   else
      printf("NO VALID USER-SUBMITTED K-CLEANUP KERNELS\n\n");
   free(dl);
   WriteMMFileWithPath(pre, "res", "UMMKCLEAN.sum", gmmb);
   return(gmmb);
}

ATL_mmnode_t *FindMUNU(ATL_mmnode_t *mb, int mu, int nu)
{
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
      if (mp->mu == mu && mp->nu == nu)
         return(mp);
   return(NULL);
}

ATL_mmnode_t *KCleanByNB
(
   int verb,
   char pre,
   ATL_mmnode_t *mmb, /* final kernels giving final supported NBs */
   ATL_mmnode_t *mkb  /* All necessary routs ku=1 to clean all kerns in mmb */
)
/*
 * Replicates mkb so that it includes all NBs in mmb, times K-clean,
 * **FREES** mkb, and returns by-NB list
 *
 * OUTPUT:
 *   <pre>UMMKCLEANBYNB.sum: non-unique K-clean for each NB in mmb
 *      mflop[1] contains estimated time for 1 K-it using K=MAX(kb/4,4)
 */
{
   ATL_mmnode_t *nkb=NULL, *mp, *np;
   int kb;
   double mf;

   printf("TIMING K-CLEAN FOR ALL SUPPORTED NBs:\n");
   for (mp=mmb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      int mb = mp->mbB, nb = mp->nbB, kb;

      kb = mp->kbB >> 2;
      kb = (kb >= 4) ? kb : 4;
      if (mp->kmaj > 1)
         kb = ((kb+mp->kmaj-1)/mp->kmaj)*mp->kmaj;

      p = FindMUNU(mkb, mp->mu, mp->nu);
/*
 *    If no user cleanup exists, generate one
 */
      if (!p)
      {
         if (mp->kmaj > 1)
            p = GetNewGenNode(pre, 0, 0, mp->mu, mp->nu, mp->kmaj, mp->kmaj);
         else
            p = GetNewGenNode(pre, 0, 0, mp->mu, mp->nu, 1, 0);
      }
      else
      {
         p = CloneMMNode(p);
         p->next = NULL;
      }
      p->nbB = nb; p->mbB = mb;  p->kbB = kb;
      p->mflop[0] = TimeMMKernel(verb, 0, p, pre, mb, nb, kb, 0, 0, mb,
                                 0, 0, -1);
      mf = mb*nb;
      p->mflop[1] = mf / p->mflop[0];   /* time in microseconds for 1 k it */
      printf("   mb=%d, nb=%d, kb=%d, mu=%d, nu=%d, mf=%.2f (%e Usec/Kit)\n",
             mb, nb, kb, p->mu, p->nu, p->mflop[0], mf);
      if (nkb)
      {
         np->next = p;
         np = p;
      }
      else
         nkb = np = p;
   }
   printf("DONE.\n\n");
   KillAllMMNodes(mkb);
   WriteMMFileWithPath(pre, "res", "UMMKCLEANBYNB.sum", nkb);
   return(nkb);
}

void TimeKClean(int verb, char pre, ATL_mmnode_t *mp)
/*
 *   mp is the ku=1, KRUNTIME K-cleanup kernel for a support NB.
 *   This routine creates an output file for supported the NB, where we
 *   document the performance for all NB different KB values.  These
 *   timings can therefore precisely document how expensive K-cleanup
 *   will be for each NB.
 *   OUTPUT:
 *   <pre>UMMKCLEAN_<nb>.TIM: timing of K-clean for nb=<nb>; there are
 *   i=nb-1 timings, mflop[0] contains time to do NB-i K its.  Will use
 *   these times to get completely accurate estimate of total time for
 *   large problems (use estimated time in CLBYNB for small probs).
 */
{
   ATL_mmnode_t *mmb, *p, *np;
   char fn[32];
   int mb = mp->mbB, nb = mp->nbB, i;

   sprintf(fn, "UMMKCLEAN_%d.TIM", mp->nbB);
   mmb = ReadMMFileWithPath(pre, "res", fn);
   if (mmb)
   {
      printf("READING IN K-CLEANUP TIMINGS FOR NB=%d:\n", nb);
      FillInGenStrings(pre, mmb);
      for (p=mmb; p; p = p->next)
      {
         int kb = p->kbB;
         assert(nb == p->nbB && mb == p->mbB);
         if (p->mflop < 0)
            p->mflop[0] = TimeMMKernel(verb, 0, p, pre, mb, nb, kb, 0, 0, kb,
                                       0, 0, -1);
         printf("   MB=%d, NB=%d, KB=%d, mu=%d, nu=%d, MFLOP=%.2f\n",
                mb, nb, kb, p->mu, p->nu, p->mflop[0]);
      }
      printf("Done.\n\n");
      WriteMMFileWithPath(pre, "res", fn, mmb);
      KillAllMMNodes(mmb);
      return;
   }

   printf("TIMING K-CLEANUP FOR MB=%d, NB=%d:\n", mb, nb);
   for (i=1; i <= nb; i++)  /* create queue of ascending KB */
   {
      p = CloneMMNode(mp);
      p->next = NULL;
      p->kbB = i;
      p->mflop[0] = TimeMMKernel(verb, 0, p, pre, mb, nb, i, 0, 0, i, 0, 0, -1);
      printf("   MB=%d, NB=%d, KB=%d, mu=%d, nu=%d, MFLOP=%.2f\n", mb, nb, i,
             p->mu, p->nu, p->mflop[0]);
      if (mmb)
      {
         np->next = p;
         np = p;
      }
      else
         mmb = np = p;
   }
   WriteMMFileWithPath(pre, "res", fn, mmb);
   KillAllMMNodes(mmb);
   printf("Done.\n");
}

void ComputeKClean(int verb, char pre, ATL_mmnode_t *mmb)
/*
 * This kernel generates K-cleanup for all routines present in mmb.
 * To do this, it finds all unique (mu,nu,kmaj) triplets, and then finds the
 * best (ku=1 || kmaj>1), K-runtime kernel for that pair (timing both user
 * and generated).
 * OUTPUT: <pre>UMMKCLEAN.sum: all unique kerns to be compiled
 *   <pre>UMMKCLEANBYNB.sum: non-unique K-clean for each NB in mmb
 *      mflop[1] contains estimated time for 1 K-it using K=MAX(kb/4,4)
 *   <pre>UMMKCLEAN_<nb>.TIM: timing of K-clean for nb=<nb>; there are
 *   i=nb-1 timings, mflop[0] contains time to do NB-i K its.  Will use
 *   these times to get completely accurate estimate of total time for
 *   large problems (use estimated time in CLBYNB for small probs).
 *
 * NOTE: we will gen/time K-cleanup kernels only in BETA=0 case, and peel
 *    the first K-block rather than the last.  This will minimize the C cost,
 *    which is more appreciable for short-K.
 */
{
   ATL_mmnode_t *mkb, *mp;
   mkb = GetUniqueKClean(verb, pre, mmb);
   mkb = KCleanByNB(verb, pre, mmb, mkb);
/*
 * For now, skip the timing until we have things working and know better
 * what timings we'll need
 */
   #if 0
/*
 * Now, create NB detailed K-cleanup files
 */
   for (mp=mkb; mp; mp = mp->next)
      TimeKClean(verb, pre, mp);
   printf("\n");
   #endif
   KillAllMMNodes(mkb);
}

void FindBestKU1
(
   int verb,
   char pre,
   int K       /* K dim, should be small, probably like 23 or 17 */
)
/*
 * Find the best possible kernel for use in low-rank update;  We only consider
 * kernels with runtime-K that handle all possible K (ku=1).  We will try
 * all legal blocking factors between 16 & 480 for this kernel, and choose
 * the one that performs best.  This kernel always used for any K not covered
 * by optimized kernels given in eAMMRES kbBs.  When we match a kbB, we
 * compare the perf of this kernel at its optimal nbB/mbB wt that of the
 * specialized kernel, and choose the best.
 *
 * OUTPUT: This routine outputs two files:
 * (1) AMMRANKK: best ku=1 kern wt best MB/NB, K=K
 * (2) AMMRANKKT: timing of this kern wt M=mbB, N=nbB, all K between 1 & maxNB
 */
{
}

double CacheRatio_all3(size_t CS, size_t mb, size_t nb, size_t kb,
                       size_t mu, size_t nu)
{ /* RETURNS: ratio of utilized cache to keep all 3 mm ops in CS */
   double dret = 1.0*kb*(mb+nb)+mb*nb;
   return(dret/CS);
}

double CacheRatio_one(size_t CS, size_t mb, size_t nb, size_t kb,
                      size_t mu, size_t nu)
{ /* RETURNS: ratio of util cache for B + working set of A/C */
   double dret = kb*nb + 2.0*(mu*kb + mu*nu);
   return(dret/CS);
}

double CacheRatio_ws(size_t CS, size_t mb, size_t nb, size_t kb,
                     size_t mu, size_t nu)
{ /* RETURNS: ratio of working set of all matmul ops to CS */
   double dret = 2.0*(mu*nu + nu*kb) + mu*kb;
   return(dret/CS);
}

typedef void (*BudgetFunc_t)(double, size_t, size_t, size_t, size_t,
                             size_t*, size_t*, size_t*);

#define MAXNB 512
void GetBlkFromBudget_all3(double thresh, size_t CS,
                           size_t mu, size_t nu, size_t ku,
                           size_t *MB, size_t *NB, size_t *KB)
{
   size_t mb=mu, nb=nu, kb=ku;
   int MGROW, NGROW, KGROW;
   do
   {
      size_t mn=mb+mu, nn=nb+nu, kn=kb+ku;
      MGROW = mn < MAXNB && (CacheRatio_all3(CS, mn, nb, kb, mu, nu) <= thresh);
      NGROW = nn < MAXNB && (CacheRatio_all3(CS, mb, nn, kb, mu, nu) <= thresh);
      KGROW = kn < MAXNB && (CacheRatio_all3(CS, mb, nb, kn, mu, nu) <= thresh);
      if (KGROW && ((!MGROW && !NGROW) || (kn <= nn && kn <= mn)))
         kb = kn;
      else if (MGROW && (!NGROW || mb < nb))
         mb = mn;
      else if (NGROW)
         nb = nn;
   }
   while (MGROW | NGROW | KGROW);
   *MB = mb;
   *NB = nb;
   *KB = kb;
}

void GetBlkFromBudget_one(double thresh, size_t CS,
                           size_t mu, size_t nu, size_t ku,
                           size_t *MB, size_t *NB, size_t *KB)
{
   size_t mb=mu, nb=nu, kb=ku;
   int NGROW, KGROW;
   do
   {
      size_t nn=nb+nu, mn=(nn/mu)*mu, mn1 = ((nn+mu-1)/mu)*mu, kn=kb+ku;
      if (mn1 - nn <= nn - mn || !mn)
         mn = mn1;
      NGROW = nn < MAXNB && (CacheRatio_one(CS, mn, nn, kb, mu, nu) <= thresh);
      KGROW = kn < MAXNB && (CacheRatio_one(CS, mb, nb, kn, mu, nu) <= thresh);
      if (NGROW && ((nn < kn && mn < kn) || !KGROW))
      {
         nb = nn;
         mb = mn;
      }
      else if (KGROW)
         kb = kn;
   }
   while (NGROW | KGROW);
   *MB = mb;
   *NB = nb;
   *KB = kb;
}

void GetBlkFromBudget_ws(double thresh, size_t CS,
                          size_t mu, size_t nu, size_t ku,
                          size_t *MB, size_t *NB, size_t *KB)
{
   size_t mb=mu, nb=nu, kb=ku;
   int KGROW;
   do
   {
      size_t kn=kb+ku, mn=(kn > mu)?(kn/mu)*mu:mu, nn=(kn>nu)?(kn/nu)*nu:nu;
      KGROW = kn < MAXNB && (CacheRatio_ws(CS, mn, nn, kn, mu, nu) <= thresh);
      if (KGROW)
      {
         mb = mn;
         nb = mn;
         kb = kn;
      }
   }
   while (KGROW);
   *MB = mb;
   *NB = nb;
   *KB = kb;
}

ATL_mmnode_t *FindBestCacheBudgetCase
(
   int verb,
   char pre,
   BudgetFunc_t GetBlocking,     /* func ptr to budget function */
   double thresh,                /* max ratio of cache to fill */
   size_t CS,                    /* size of cache we are optimizing for */
   int imf,                      /* entry in mflop[] to use */
   ATL_mmnode_t *mmb             /* list of cases to try */
)
/*
 * RETURNS: clone of best-peforming kernel in mmb for kb=kb, mb & nb
 *          near-square and within budget
 */
{
   ATL_mmnode_t *mmB=NULL, *mp, *p;
   double mf, mfB=0.0;

   printf("Finding best case for cache budget case=%d, CS=%.0f elts\n",
          imf, CS*thresh);
   for (mp=mmb; mp; mp = mp->next)
   {
      size_t mb, nb, kb;
      GetBlocking(thresh, CS, mp->mu, mp->nu, mp->ku, &mb, &nb, &kb);
      p = CloneMMNode(mp);  /* can't use mp, since may switch KRUNTIME */
      mf = TimeMMKernel_KB(verb, 0, p, pre, mb, nb, kb, 1, 0, -1);
      printf("   ID=%d, mb=%d, nb=%d, kb=%d, RTK=%d, MFLOP=%.2f\n", p->ID,
             (int)mb, (int)nb, (int)kb, FLAG_IS_SET(p->flag, MMF_KRUNTIME), mf);
      if (mf > mfB)
      {
         if (mmB)
            KillMMNode(mmB);
         p->mbB = mb;
         p->nbB = nb;
         p->kbB = kb;
         mmB = p;
         mfB = mmB->mflop[imf] = mf;
      }
      else
         KillMMNode(p);
   }
   printf("BEST CASE %s: mb=%d, nb=%d, kb=%d, RTK=%d, MFLOP=%.2f\n\n",
          mmB->rout ? mmB->rout : "GENNED",
          mmB->mbB, mmB->nbB, mmB->kbB, FLAG_IS_SET(mmB->flag, MMF_KRUNTIME),
          mmB->mflop[imf]);
   mmB->next = NULL;
   return(mmB);
}

ATL_mmnode_t *FindBestCacheBudgetCases
(
   int verb,
   char pre,
   size_t CS,                    /* size of cache we are optimizing for */
   ATL_mmnode_t *mmb             /* list of cases to try */
)
/*
 * This case attempts to find the best kernel for 3 cases of interest:
 * (1) All 3 matrices fit in CS -- this case is designed for when we wish
 *     to reuse at least one of the matrices *across* mmkern calls.  It is
 *     particularly good for complex arithmetic, or when CS is large enough
 *     that A&B are reused so much internally to a mmkern call that it makes
 *     sense to retain C in cache for the next mmkern call in K-loop.
 * (2) All of B fits in CS, and so does the working set of A&C.  This one
 *     reuses all internal ops from L1, but won't allow any full op reuse
 *     across multiple GEMM calls.  Usually best for small-to-medium cache
 *     sizes.
 * (3) All of working set of A/B/C fit in cache.  This case provides maximal
 *     NB, where only the mu*KB panel of A is reused from the L1 internally
 *     to the algorthm.  It is essentially an L2-blocked algorithm internally,
 *     but can be useful on those archs where the best sustained bandwidth
 *     comes from one L1 load (A) and 1 L2 load (B).
 */
{
   ATL_mmnode_t *mm3, *mm1, *mmw;
   mm3 = FindBestCacheBudgetCase(verb, pre, GetBlkFromBudget_all3, 1.0,
                                 CS, 1, mmb);
   mm1 = FindBestCacheBudgetCase(verb, pre, GetBlkFromBudget_one, 1.0,
                                 CS, 2, mmb);
   mmw = FindBestCacheBudgetCase(verb, pre, GetBlkFromBudget_ws, .90,
                                 CS, 3, mmb);
   mm3->next = mm1;
   mm1->next = mmw;
   mmw->next = NULL;

   printf("3CASES: ID=%d/%d/%d, mb=%d/%d/%d, nb=%d/%d/%d, kb=%d/%d/%d\n",
          mm3->ID, mm1->ID, mmw->ID, mm3->mbB, mm1->mbB, mmw->mbB,
          mm3->nbB, mm1->nbB, mmw->nbB, mm3->kbB, mm1->kbB, mmw->kbB);
   printf("        RTK=%d/%d/%d, MFLOP=%.2f/%2f/%2f\n\n",
          FLAG_IS_SET(mm3->flag, MMF_KRUNTIME),
          FLAG_IS_SET(mm1->flag, MMF_KRUNTIME),
          FLAG_IS_SET(mmw->flag, MMF_KRUNTIME),
          mm3->mflop[1], mm1->mflop[2], mmw->mflop[3]);

   return(mm3);
}

int KernelIsSame(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if kernels are the same except for blocking, 0 otherwise
 */
{
/*
 * Two generated kernels are the same if mu,nu,ku,VLEN,flag are the same.
 * NOTE: if we make generator handle muladd, etc, MUST UPDATE HERE!!!
 */
   if (p0->ID == 0 && p1->ID == 0)
      return(p0->mu == p1->mu && p0->nu == p1->nu && p0->ku == p1->ku &&
             p0->vlen == p1->vlen && p0->flag == p1->flag &&
             p0->kmaj == p1->kmaj);
/*
 * If both are user kernels, then they may be repeats.  For user kernels,
 * they are the same if both ID and flag match, else they are not.
 */
   else if (p0->ID > 0 && p1->ID > 0)
      return(p0->ID == p1->ID && p0->flag == p1->flag);
   return(0);  /* Can't be the same if above criteria fails */
}

int KernelIsUnique(ATL_mmnode_t *mmb, ATL_mmnode_t *mmp)
/*
 * Determines if mmp is the first mention of a unique kernel in mmb, or not.
 * For user cases (ID > 0), (ID,flag) together make a unique kernel.
 * For user generated cases, if they match on : mu,nu,ku,VLEN,flag
 *
 * RETURNS: 0 if mmp appears in mmb before mmp, else 1
 */
{
   ATL_mmnode_t *mp;
   if (mmp == mmb)
      return(1);
   for (mp=mmb; mp && mp != mmp; mp = mp->next)
      if (KernelIsSame(mmp, mp))
         return(0);
   return(1);  /* didn't find it, must be first time in list */
}
/*
 * Returns a non-repetitive list of user kernels (ID>0) found in rb.  Note that
 * differing compilations of the same kernel are reduced to one entry.
 * rb is left unchanged.
 */
ATL_mmnode_t *GetUniqueUserKerns(ATL_mmnode_t *rb)
{
   ATL_mmnode_t *ub=NULL, *p;

   if (!rb)
      return(NULL);
   for (p=rb; p; p = p->next)
      if (p->ID > 0)
         break;
   if (!p)
      return(NULL);
   ub = CloneMMNode(p);
   for (p=p->next; p; p = p->next)
   {
       if (p->ID > 0)
       {
          ATL_mmnode_t *np;
          int ID = p->ID;

          for (np=ub; np; np = np->next)
             if (np->ID == ID)
                break;
          if (!np)
          {
             np = CloneMMNode(p);
             np->next = ub;
             ub = np;
          }
       }
   }
   return(ub);
}

ATL_mmnode_t *TimeKBRegion
(
   int verb,
   char pre,
   ATL_mmnode_t *mmk,            /* kernel to time throughout region */
   int kbmin,                    /* start of region */
   int kbend,                    /* largest kb in region */
   int kincD                     /* default stride between kernel timings */
)
/*
 * Returns list of timings of kernel mmk using near-square cases with KB
 * varying between kbmin - kbend.  All cases that are legal and incremented
 * by kinc are tried, as are all perfectly square cases
 */
{
   ATL_mmnode_t *mmb=NULL, *mp, *mpB=NULL;
   const int ku = mmk->ku, mu=mmk->mu, nu=mmk->nu;
   int kstart, kinc, kend, k, ksq, ksqinc;
   double mf, mfB=0.0;
/*
 * Get starting and ending point that is legal for this kernel.
 */
   kstart = Mmax(mmk->kbmin, kbmin);
   kstart = ((kstart+ku-1)/ku)*ku;
   kend = ((kbend+ku-1)/ku)*ku;
   if (mmk->kbmax)
      kend = Mmin(kend, mmk->kbmax);
   k = kstart;
/*
 * square inc always lcm(mu,nu,ku).  Normal increment is always at least
 * as big as the default stride, but must be a multiple of the kernel's ku
 */
   ksqinc = Mylcm(mu, nu);
   ksqinc = Mylcm(ksqinc, ku);
   for (kinc=ku; kinc < kincD; kinc += ku);
   if (kstart <= kend)
   {
      int kb = k;
      printf("TIMING %s mu=%d, mu=%d, For KB=[%d,%d]:\n",
             mmk->rout ? mmk->rout : "Genkern", mu, nu, kstart, kend);
      ksq = ((kstart+ksqinc-1)/ksqinc)*ksqinc;
      do
      {
         ATL_mmnode_t *p;
         const int mb=((kb+mu-1)/mu)*mu, nb=((kb+nu-1)/nu)*nu;

         p = CloneMMNode(mmk);
         p->mbB = mb; p->nbB = nb; p->kbB = kb;
         mf = TimeMMKernel_KB(verb, 0, p, pre, mb, nb, kb, 1, 0, -1);
         printf("   mb=%d, nb=%d, kb=%d, KRUN=%d, MFLOP=%.2f\n",
                mb, nb, kb, FLAG_IS_SET(p->flag, MMF_KRUNTIME), mf);
         p->mflop[0] = mf;
         if (mf > mfB)
         {
            mfB = mf;
            mpB = p;
         }
         if (mmb)
         {
            mp->next = p;
            mp = p;
         }
         else
            mmb = mp = p;
         if (kb == ksq)
            ksq += ksqinc;
         if (kb == k)
            k += kinc;
         kb = Mmin(k,ksq);
      }
      while(kb <= kend);
      printf("DONE, best case %s mb=%d, nb=%d, kb=%d, MFLOP=%.2f\n\n",
             mmk->rout ? mmk->rout : "Genkern",
             mpB->mbB, mpB->nbB, mpB->kbB, mfB);
   }
   else
   {
      printf("KERNEL %s mu=%d, mu=%d, has no legal cases in KB=[%d,%d]!\n\n",
             mmk->rout ? mmk->rout : "Genkern", mu, nu, kstart, kend);
   }
   return(mmb);
}

ATL_mmnode_t *TimeAllKBRegions
(
   int verb,
   char pre,
   ATL_mmnode_t *mmk,            /* kernel to time throughout regions */
   int kb1,                      /* rough kb ending region 1 */
   int kb2,                      /* rough kb ending region 2 */
   int kb3                       /* maxKB to ever try */
)
/*
 * Times mmk for KBs in the three regions, returns list
 */
{
   const int ku = mmk->ku;
   ATL_mmnode_t *mmb, *mp;

   mmb = TimeKBRegion(verb, pre, mmk, 24, kb3, 4);
   return(mmb);
}

ATL_mmnode_t *FindCacheBudgetCasesByKB
(
   int verb,
   char pre,
   size_t CS,                    /* size of cache we are optimizing for */
   ATL_mmnode_t *mmb             /* list of cases to try */
)
/*
 * This routine is responsible for:
 * (1) Find the best performing kernels out of mmb for our 3 cache budget
 *     cases
 * (2) Free mmb
 * (3) For each unique kernel, find perf of kernel for all supported KBs
 *     in the budgetary regions
 * (4) Merge these lists, and winnow underperforming cases
 * (5) RETURN: queue of all supported KBs
 */
{
   ATL_mmnode_t *mm3b, *mp;
/*
 * See if we just need to rerun cases
 */
   mm3b = ReadMMFileWithPath(pre, "res", "bAMMRES.sum");
   if (mm3b)
   {
      int i=0;
      printf("READING IN LARGE KERNEL CASES FROM res/<pre>bAMMRES:\n");
      FillInGenStrings(pre, mm3b);
      for (mp=mm3b; mp; mp = mp->next)
      {
         if (mp->mflop[0] <= 0.0)
            mp->mflop[0] = TimeMMKernel(verb, 0, mp, pre, mp->mbB, mp->nbB,
                                        mp->kbB, 0, 0, 0, 1, 0, -1);
         printf("   ID=%d, %s: MB=%d, NB=%d, KB=%d, KRUN=%d, MFLOP=%.2f\n",
                mp->ID, mp->rout ? mp->rout : "Gennedkern",
                mp->mbB, mp->nbB, mp->kbB, FLAG_IS_SET(mp->flag, MMF_KRUNTIME),
                mp->mflop[0]);
         i++;
      }
      printf("DONE %d CASES.\n\n", i);
      return(mm3b);
   }
/*
 * Find best performing kernels for each of our 3 cache budgets
 */
   mm3b = FindBestCacheBudgetCases(verb, pre, CS, mmb);
   KillAllMMNodes(mmb);
/*
 * Get list of performance of all-3 in-cache kernel in all 3 cache regions
 */
   mmb = TimeAllKBRegions(verb, pre, mm3b, mm3b->kbB, mm3b->next->kbB,
                          mm3b->next->next->kbB);
   for (mp=mm3b->next; mp; mp = mp->next)
   {
      if (KernelIsUnique(mm3b, mp))
      {
         ATL_mmnode_t *p, *p2;
         p = TimeAllKBRegions(verb, pre, mp, mm3b->kbB, mm3b->next->kbB,
                              mm3b->next->next->kbB);
         p2 = MergeCases(0, mmb, p);
         KillAllMMNodes(p);
         KillAllMMNodes(mmb);
         mmb = p2;
      }
   }
   KillAllMMNodes(mm3b);
/*
 * Now, get rid of any blocking factor that is slower than the preceeding one
 */
   mmb = WinnowCases(0, mmb);
   WriteMMFileWithPath(pre, "res", "bAMMRES.sum", mmb);
   return(mmb);
}

ATL_mmnode_t *DecentGenCase(int verb, char pre, int nreg)
{
/*
 * Find out what vectorization, if any, to use in generating kernels
 */
   int mu, nu;
   SetGenVec(verb, pre);
/*
 * 120 = LCM(2,3,4,5,6,8), and large enough to stress mu/nu
 */
   return(FindDefMUNU(verb, pre, nreg, 0, 120, 1, &mu, &nu));
}

/*
 * This routine finds kernels to use in low-rank-K update. For 3 <= K <= 15,
 * it tries all kernels and chooses the best performing; In this search
 * we consider only compile-time K kernels, since runtime kernels will be
 * selected by general (K>15) search.
 * (K=1 and K=2 are handled by GER and GER2).
 * We consider only user-generated kernels; for any problem sizes that are
 * not supported, we will use the normal K or K-clean routines.  This list
 * is just to allow for hand-tuning small-K special cases.
 * This routine produces output file <pre>AMMLOWK.sum
 */
ATL_mmnode_t *GetLowRankKKernels
(
   int verb,            /* verbosity */
   char pre,            /* s,d */
   int MB,              /* default mb to time with */
   int NB,              /* default nb to time with */
   ATL_mmnode_t *inb    /* all working ukerns */
)
{
   int k, ik;
   ATL_mmnode_t *rkKb=NULL, *mp;
/*
 * Get rid of all K-runtime kernels from consideration
 */
   while (inb && FLAG_IS_SET(inb->flag, MMF_KRUNTIME))
      inb = KillMMNode(inb);
   if (inb)
   {
      ATL_mmnode_t *prev=inb;
      mp = inb->next;
      while (mp)
      {
         if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
         {
            mp = KillMMNode(mp);
            prev->next = mp;
         }
         else
         {
            prev = mp;
            mp = mp->next;
         }
      }
   }
   else
      return(NULL);
   for (ik=0; ik < 2; ik++)
   {
      int kbeg, kend, kinc;
      if (!ik)
      {
         kbeg = 96;
         kend = 16;
         kinc = 16;
      }
      else
      {
         kbeg = 15;
         kend = 3;
         kinc = 1;
      }
      for (k=kbeg; k >= kend; k -= kinc)
      {
         printf("FINDING BEST USER-PROVIDED KERNEL FOR K=%d:\n", k);
         ATL_mmnode_t *best=NULL;
         for (mp=inb; mp; mp = mp->next)
         {
            const int mu = mp->mu, nu = mp->nu, ku = mp->ku;
            const int mb = (MB/mu)*mu, nb = (NB/nu)*nu;
            const int KK = (mp->kmaj < 2) ? k : ((k+ku-1)/ku)*ku;
            double mf;

            assert(mb && nb);
            if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME) || KK%ku)
            {
               printf("   skipping %d. %s, KRUN=%d, ku=%d\n", mp->ID, mp->rout,
                      FLAG_IS_SET(mp->flag, MMF_KRUNTIME), ku);
               continue;
            }
            if ((mp->kbmin && k < mp->kbmin) || (mp->kbmax && KK > mp->kbmax))
            {
               printf("   skipping %d. %s, kbmin,max=%d,%d, K=%d\n",
                      mp->ID, mp->rout, mp->kbmin, mp->kbmax, KK);
               continue;
            }
            mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, KK, 0, 0, 0, 0, 0, -1);
            if (KK != k)
               mf = (mf*k) / (double)KK;
            printf("   %d. %s: mb=%d, nb=%d, MFLOP=%.2f\n", mp->ID, mp->rout,
                   mb, nb, mf);
            if (!best)
            {
               best = mp;
               mp->mflop[0] = mf;
               mp->mbB = mb;
               mp->nbB = nb;
               mp->kbB = k;
            }
            else if (best->mflop[0] < mf)
            {
               best = mp;
               mp->mflop[0] = mf;
               mp->mbB = mb;
               mp->nbB = nb;
               mp->kbB = k;
            }
         }
         if (best)
         {
            best = CloneMMNode(best);
            best->next = rkKb;
            rkKb = best;
            printf("BEST FIXED-%d KERNEL: %d. %s MFLOP=%.2f\n\n",
                   k, best->ID, best->rout, best->mflop[0]);
         }
         else
            printf("NO SPECIAL CASE for K=%d\n\n", k);
      }
   }
   KillAllMMNodes(inb);
   return(rkKb);
}

/*
 * Finds list of best run-time kernel, ranked by KU.  Higher KUs are not
 * retained unless they beat any lower-KU kernel that divides that KU evenly
 */
ATL_mmnode_t *GetRuntimeKKernels
(
   int verb,            /* verbosity */
   char pre,            /* s,d */
   int MB,              /* default mb to time with */
   int NB,              /* default nb to time with */
   ATL_mmnode_t *inb    /* all working ukerns */
)
{
   ATL_mmnode_t *mp;
   int KU;
/*
 * Get rid of all non-K-runtime kernels from consideration
 */
   while (inb && !FLAG_IS_SET(inb->flag, MMF_KRUNTIME))
      inb = KillMMNode(inb);
/*
 * Got rid of any at base, now get rid of non-K-runtime from internal nodes
 */
   if (inb)
   {
      ATL_mmnode_t *prev=inb;
      mp = inb->next;
      while (mp)
      {
         if (!FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
         {
            mp = KillMMNode(mp);
            prev->next = mp;
         }
         else
         {
            prev = mp;
            mp = mp->next;
         }
      }
   }
   else
      return(NULL);
   KU = inb->ku;
   for (mp=inb->next; mp; mp = mp->next)
      KU = Mylcm(KU, mp->ku);
   if (KU > 32)
      KU = 32;
   else
      KU = ((16+KU-1)/KU)*KU;
   printf("TRYING ALL RUNTIMEK KERNS WITH MB=%d, NB=%d, KB=%d:\n", MB, NB, KU);
   for (mp=inb; mp; mp = mp->next)
   {
      int mu=mp->mu, nu=mp->nu, ku=mp->ku;
      int mb = (MB/mu)*mu, nb = (NB/nu)*nu, kb = (KU/ku)*ku;
      double mf;
      assert(mb && nb && kb);
      mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, kb, 0, 0, 0, 0, 0, -1);
      printf("   %d. %s: mb=%d, nb=%d, MFLOP=%.2f\n", mp->ID, mp->rout,
             mb, nb, mf);
      mp->mflop[0] = mf;
      mp->mbB = mb;
      mp->nbB = nb;
      mp->kbB = kb;
   }
   printf("\n");
   inb = ATL_SortMMNodesByMflop(0, inb);
   if (inb->ku == 1)
   {
      KillAllMMNodes(inb->next);
      inb->next = NULL;
   }
   else
   {
      ATL_mmnode_t *p;
/*
 *    Go thru sorted list, and kill all slower nodes that don't add new K
 */
      for (p=inb; p; p = p->next)
      {
         ATL_mmnode_t *prev=p;
         mp = p->next;
         while (mp)
         {
            if (mp->ku % p->ku == 0)
            {
               mp = KillMMNode(mp);
               prev->next = mp;
            }
            else
            {
               prev = mp;
               mp = mp->next;
            }
         }
      }
   }
   if (!inb)
      printf("NO RETAINED RUNTIME KERNELS.\n\n");
   else
   {
      printf("RETAINED RUNTIME KERNELS:\n");
      for (mp=inb; mp; mp = mp->next)
         printf("   %d. %s: ku=%d, MFLOP=%.2f\n", mp->ID, mp->rout, mp->ku,
                mp->mflop[0]);
      printf("DONE.\n");
   }
   return(inb);
}
/*
 * RETURNS: 1 if mmc is slower than any kernel in mmb
 */
int IsSlowerThanList
(
   int verb,            /* verbosity */
   char pre,            /* s,d */
   int MB,
   int NB,              /* default mb/nb to time with */
   ATL_mmnode_t *mmc,  /* candidate mmkern */
   ATL_mmnode_t *mmb   /* kernels to time candidate against */
)
{
   ATL_mmnode_t *mp;
   double mfc, mf;
   int mu, nu, ku;
   int mb, nb, kb, KB;

   if (!mmb)
      return(0);
   kb = mmc->kbB;
   mu = mmc->mu;
   nu = mmc->nu;
   mb = (MB/mu)*mu;
   nb = (NB/nu)*nu;
   assert(mb && nb && kb);
   KB = (mmc->kmaj < 2) ? kb : ((kb+mmc->ku-1)/mmc->ku)*mmc->ku;
   mfc = TimeMMKernel(verb, 0, mmc, pre, mb, nb, KB, 0, 0, 0, 0, 0, -1);
   mfc = (kb*mfc)/(double)KB;
   mmc->mflop[1] = mfc;
   kb = mmc->kbB;
   for (mp=mmb; mp; mp = mp->next)
   {
      ku = mp->ku;
      if (mp->kbmin && kb < mp->kbmin)
         continue;
      if (mp->kbmax && kb > mp->kbmax)
         continue;
      if (kb%ku == 0 || mp->kmaj > 1)
      {
         int KK = (mp->kmaj > 1) ? ((kb+ku-1)/ku)*ku : kb;
         mu = mp->mu;
         nu = mp->nu;
         mb = (MB/mu)*mu;
         nb = (NB/nu)*nu;
         assert(mb && nb);
         mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, KK, 0, 0, 0, 0, 0, -1);
         mp->mflop[1] = (kb*mfc)/(double)KK;
         if (mf > mfc)
         {
            return(1);
         }
      }
   }
   return(0);
}

/*
 * Finds best-performing square cases in list of pre-existing kernels, mmb.
 * Does not modify original mmb list, and will return only kernels that are
 * faster than smaller square cases.
 * RETURNS: new list of all square cases that got best performance from
 *          original mmb.
 */
ATL_mmnode_t *FindBestSquareCases(char pre, int verb, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mmp, *mmSQ=NULL, *prev=NULL;
   int maxNB=0, maxU=0, i;

   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      if (mmp->mu > maxU)
         maxU = mmp->mu;
      if (mmp->nu > maxU)
         maxU = mmp->nu;
      if (mmp->nbB > maxNB)
         maxNB = mmp->nbB;
      if (mmp->mbB > maxNB)
         maxNB = mmp->mbB;
      if (mmp->kbB > maxNB)
         maxNB = mmp->kbB;
   }
   maxNB = ((maxNB+maxU-1) / maxU)*maxU;

   for (i=4; i < maxNB; i++)
   {
      mmp = BestForThisNB(verb, pre, mmb, i, i-1, i+1, 1);
      if (mmSQ)
      {
         if (prev->mflop[0] >= mmp->mflop[0])
            KillMMNode(mmp);
         else
         {
            prev->next = mmp;
            prev = mmp;
         }
      }
      else
         mmSQ = prev = mmp;
   }
   return(mmSQ);
}
/*
 * ASSUMES: mmb comes in ordered from smallest blocking to largest.
 * Kills any kernel that is not faster than a smaller-NB kernel.
 */
ATL_mmnode_t *HarshPrune(ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mmp=mmb;
   while (mmp)
   {
      if (!mmp->next)
         return(mmb);
      do
      {
         if (!mmp->next)
            return(mmb);
         if (mmp->mflop[0] > mmp->next->mflop[0])
            mmp->next = KillMMNode(mmp->next);
         else
            mmp = mmp->next;
      }
      while (mmp);
   }
}
/*
 * ASSUMES: mmb comes in ordered from smallest blocking to largest.
 * Kills any kernel that is not faster than a smaller-NB kernel.
 */
ATL_mmnode_t *TolerantPrune
(
   ATL_mmnode_t *mmb,   /* list of kernels to prune */
   double tol           /* set > 1 to allow slow kernels to stay */
)                       /* < 1 to make pruning harsher */
{
   ATL_mmnode_t *mmp=mmb;
   double mfB = 0.0;
   while (mmp)
   {
      if (!mmp->next)
         return(mmb);
      do
      {
         if (!mmp->next)
            return(mmb);
         if (mmp->mflop[0] > mfB)
            mfB = mmp->mflop[0];
         if (mfB > tol*mmp->next->mflop[0])
            mmp->next = KillMMNode(mmp->next);
         else
            mmp = mmp->next;
      }
      while (mmp);
   }
}
#define CON_NOKVEC 0
#define CON_NOMVEC 1
#define CON_NOCOMPK 2
/*
 * RETURNS: N-length queue wt best-performing kernel for specified dims
 */
ATL_mmnode_t *FindBestKerns
(
   char pre,
   int verb,
   ATL_mmnode_t *mmb,            /* candidate user-supplied kerns */
   int N,                        /* number of forced block factors */
   int *mbs, int *nbs, int *kbs, /* forced block factors */
   int C,                        /* constraints */
   float tol                     /* tolerance in getting rid of slow kerns */
)
{
   ATL_mmnode_t *nmmb=NULL, *mp, *mpB, *mpG, *mpg;
   const int TRYKVEC=!(C & (1<<CON_NOKVEC));
   const int TRYMVEC=!(C & (1<<CON_NOMVEC));
   int i, mbB=0;
   double mfB, mf;
/*
 * Find basic register blocking for generated case
 */
   mpG = DecentGenCase(verb, pre, 32);
   if (mpG->vlen != VLEN[VECi])
      VECi = VTSC;
   assert (mpG->kmaj < 2);
/*
 * Loop over all required blocking factors
 */
   printf("SEARCHING %d USER BLOCKINGS:\n", N);
   for (i=N-1; i >= 0; i--)
   {
      int mb=mbs[i], nb=nbs[i], kb=kbs[i], MB=mbs[i];
      int kuG, KK=kb, veci0=VECi;
      double mfK=0.0, mfM=0.0;
      int iveci = VECi, VL=VLEN[VECi];
/*
 *    If mb not specified, choose one near specified KB
 */
      printf("\n   FINDING KERNEL FOR MB=%d, NB=%d, KB=%d:\n", mb, nb, kb);
/*
 *    Find best of M- or K-vectorized generated code for this problem size
 */
      mpg = mpB = NULL;
      if (TRYMVEC)
      {
         int muG=mpG->mu, nuG=mpG->nu, ut=muG*nuG;
/*
 *       If mb is not specified, make it a multiple of mu & vlen
 */
         if (!MB)
         {
            int mb0;
            mb = Mylcm(muG, VL);
            mb0 = (kb/mb)*mb;
            mb = ((kb+mb-1)/mb)*mb;
            if (mb0 && mb-kb > kb-mb0)
               mb = mb0;
         }
/*
 *       If specified mb not a multiple of vlen, reduce vlen until it is
 */
         else if (mb%VL)
         {
            if (VECi == VTAVX && !(mb%VLEN[VTSSE]))
               VECi = VTSSE;
            else
               VECi = VTSC;
         }
/*
 *       Make register block a multiple of mandated block
 */
         muG = (muG / VL)*VLEN[VECi];
         ut = muG * nuG;
         while (muG > 1 && mb%muG)
            muG -= VLEN[VECi];
         if (muG < 1)
            muG = 1;
         while(muG*(nuG+1) <= ut)
            nuG++;
         while (nb%nuG)
            nuG--;
         kuG = 1;
         mpg = GetNewGenNode(pre, kb, 0, muG/VLEN[VECi], nuG, 1, 0);
         mpg->mbB = mb;
         mfM = TimeMMKernel(verb, 0, mpg, pre, mb, nb, KK, KK, KK, mb, 1, 0,-1);
         printf("      Gen mu=%d, nu=%d, kmaj=%d, mf=%.2f\n", muG, nuG, 0, mfM);
      }
      VECi = iveci;
/*
 *    Generator currently does not support k-vectorized kernels
 */
      if (0 && TRYKVEC)
      {
         int muG, nuG=mpG->nu, ut;

         if (VLEN[VECi] > 1)
            muG = mpG->mu / VLEN[VECi];
         else
            muG = mpG->mu;
         ut = muG * nuG;
/*
 *       If kb not a multiple of vlen, reduce it until it is
 */
         if (kb%VLEN[VECi])
         {
            if (VECi == VTAVX && !(kb%VLEN[VTSSE]))
               VECi = VTSSE;
            else
               VECi = VTSC;
         }
         kuG = VLEN[VECi];
/*
 *       Make register block a multiple of mandated block
 */
         while (mb%muG)
            muG--;
         while(muG*(nuG+1) <= ut)
            nuG++;
         while (nb%nuG)
            nuG--;
         if (kuG > 1)
         {
            mpB =  GetNewGenNode(pre, kb, 0, muG, nuG, kuG, kuG);
            mpg->mbB = mb;
            mfK = TimeMMKernel(verb, 0, mpB, pre, mb, nb, KK, KK, KK, mb,
                               1, 0,-1);
            printf("      Gen mu=%d, nu=%d, kmaj=%d, mf=%.2f\n", muG, nuG, kuG,
                   mfK);
         }
      }
      if (mfM >= mfK)  /* use M-vectorized gen kernel */
      {
         if (mpB)
            KillMMNode(mpB);
         mpB = mpg;
         mfB = mfM;
      }
      else
      {
         if (mpg)
            KillMMNode(mpg);
         mpg = mpB;
         mfB = mfK;
      }
      VECi = iveci;
      if (KK != kb)
         mfB *= (double)kb / (double)KK;
      mpg->mflop[0] = mfB;
/*
 *    Search through all user kernels, and take best performing
 */
      mpB = mpg;
      mbB = mpg->mbB;
      for (mp=mmb; mp; mp = mp->next)
      {
         const int mu=mp->mu, nu=mp->nu, ku=mp->ku, kmaj=mp->kmaj;
         const int KB = (kmaj < 2) ? kb : ((kb+ku-1)/ku)*ku;
         int SKIP=0;
         if (!MB)
         {
            if (nb%nu || KB%ku)
               SKIP=1;
            else
            {
               int mb0;
               mb0 = (kb/mu)*mu;
               mb = ((kb+mu-1)/mu)*mu;
               if (mb0 && (mb-kb > kb-mb0))
                  mb = mb0;
            }
         }
         else if (mb%mu || nb%nu || KB%ku)
            SKIP=1;
         if (mp->kbmin && mp->kbmin > kb)
            SKIP=1;
         else if (mp->kbmax && mp->kbmax < kb)
            SKIP=1;
         if (SKIP)
         {
            printf("      %d-%s: SKIPPED\n", mp->ID, mp->rout);
            continue;
         }
         mf = TimeMMKernel(verb, 0, mp, pre, mb, nb, KB, KB, KB, mb, 1, 0, -1);
         if (KB != kb)
            mf *= (double)kb / (double)KB;
         printf("      %d-%s: %.2f\n", mp->ID, mp->rout, mf);
         if (mf > mfB)
         {
            mfB = mf;
            mpB = mp;
            mbB = mb;
         }
      }
      printf("    [M,N,K]B=%d,%d,%d BEST: %d-%s %.2f\n", mbB, nb, kb,
             mpB->ID, mpB->rout, mfB);
      if (mpB == mpg)
      {
         mpg->next = nmmb;
         nmmb = mpg;
      }
      else
      {
         KillMMNode(mpg);
         mp = CloneMMNode(mpB);
         mp->next = nmmb;
         nmmb = mp;
      }
      nmmb->mflop[0] = mfB;
      nmmb->mbB = mbB;
      nmmb->nbB = nb;
      nmmb->kbB = kb;
      VECi = veci0;
   }
   KillAllMMNodes(mmb);
   KillMMNode(mpG);
   if (tol > 0.0)
      nmmb = TolerantPrune(nmmb, tol);
   WriteMMFileWithPath(pre, "res", "uAMMFRC.sum", nmmb);
   printf("DONE.\n");
   return(nmmb);
}

/*
 * This routine applies the contraints C to mmb, and returns winnowed queue
 * RETURNS: Q wt all kernels failing constraints stropped out
 */
ATL_mmnode_t *ApplyConstraints(ATL_mmnode_t *mmb, int C)
{
    if (C)
    {
       ATL_mmnode_t *mp = mmb, *next;
       while (mp)
       {
          next = mp->next;
          if ((C & (1<<CON_NOKVEC)) && mp->kmaj > 1)
             mmb = KillMMNodeFromQ(mmb, mp);
          else if ((C & (1<<CON_NOMVEC)) && mp->kmaj < 2)
             mmb = KillMMNodeFromQ(mmb, mp);
          else if ((C & (1<<CON_NOCOMPK)) &&
                   !FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
             mmb = KillMMNodeFromQ(mmb, mp);
          mp = next;
       }
    }
    return(mmb);
}
void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr,
   "   -S <file> : take kerns from sum file, force square sizes\n");
   fprintf(stderr, "   -T #  : set pruning tolerance for slower kernels:\n");
   fprintf(stderr,
   "      # <= 0.0 : keep all legal kernel sizes regardless of performance\n");
   fprintf(stderr,
   "      # >=0.0 : delete all kerns wt perf*# < maxSmallerPerf\n");
   fprintf(stderr,
   "                maxSmallerPerf= max perf found in smaller-sized kern\n");
   fprintf(stderr, "   -p [s,d]: set precision prefix (d) \n");
   fprintf(stderr, "   -b # nb1 ... nb# : square NBs to force\n");
   fprintf(stderr, "   -B # mb1 nb1 kb1 ... mb# nb# kb#: dims to force\n");
   fprintf(stderr,
           "   -M # abc1 ... abc#: which matblks should be cache flushed\n");
   fprintf(stderr, "   -v <verb> : set verbosity (1)\n");
   fprintf(stderr, "   -C [kmc] : Constrain kernel choice:\n");
   fprintf(stderr, "      k : don't allow K-vectorized storage\n");
   fprintf(stderr, "      m : don't allow M-vectorized storage\n");
   fprintf(stderr, "      c : don't allow compile-time K kernels\n");
   fprintf(stderr, "   -K <1/0> : do/don't generate K cleanup\n");

   fprintf(stderr,
      "   -F <b0> <bN> <KI> <idx> <rfn>: brute-force blocking search:\n");
   fprintf(stderr, "       b0: smallest value to try for all dims\n");
   fprintf(stderr, "       b1: largest value to try for all dims\n");
   fprintf(stderr, "       KI: min K increment\n");
   fprintf(stderr, "      idx: index in rfn to use; -1 means last\n");
   fprintf(stderr, "      rfn: search result file name to read kern from\n");
   exit(ierr ? ierr : -1);
}

ATL_mmnode_t *GetFlags(int nargs, char **args, char *PRE, int *VERB, int *NN,
                       int *C, int **MBS, int **NBS, int **KBS, float *TOL,
                       int *KCLEAN)
{
   ATL_mmnode_t *mmb=NULL;
   int B0=0, BN, KI;
   int *mbs=NULL, *nbs=NULL, *kbs=NULL;
   int i, k, j, N=0, ALL=0, MV;
   char *cs;

   *TOL = 0.0;
   *VERB = 1;
   *PRE = 'd';
   *C = 0;
   *KCLEAN = 1;
   for (i=1; i < nargs; i++)
   {
      int n;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'T':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *TOL = atof(args[i]);
         if (*TOL < 0.0)
            *TOL = 0.0;
         break;
      case 'F':  /* <b0> <bN> <KI> <idx> <rfn> */
         if (i+5 >= nargs)
            PrintUsage(args[0], i-1, NULL);
         else
         {
            int I, k;
            ATL_mmnode_t *mp;

            B0 = atoi(args[i+1]);
            BN = atoi(args[i+2]);
            KI = atoi(args[i+3]);
            I = atoi(args[i+4]);
            mmb = ReadMMFile(args[i+5]);
            assert(mmb);
            if (I < 0)
               for (mp=mmb; mp->next; mp = mp->next);
            else
               for (k=0, mp=mmb; k < I && mp; k++, mp = mp->next);
            assert(mp);
            mp = CloneMMNode(mp);
            KillAllMMNodes(mmb);
            mmb = mp;
         }
         i += 5;
         break;
      case 'M':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         n = atoi(args[i]);
         for (MV=k=0; k < n; k++)
         {
           if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            if (args[i][0] == 'a' || args[i][0] == 'A')
               MV |= 1;
            else if (args[i][0] == 'b' || args[i][0] == 'B')
               MV |= 2;
            else if (args[i][0] == 'c' || args[i][0] == 'C')
               MV |= 4;
            else
               PrintUsage(args[0], -i, "UNKNOWN MATRIX FOR -M");
         }
         IMVS = MV;
         break;
      case 'S':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         mmb = ReadMMFile(args[i]);
         break;
      case 'K':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *KCLEAN = atoi(args[i]);
         break;
      case 'v':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *VERB = atoi(args[i]);
         break;
      case 'p':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *PRE = tolower(args[i][0]);
        assert(*PRE == 's' || *PRE == 'd' || *PRE == 'z' || *PRE == 'c');
        break;
      case 'C':  /* -C <constraint string */
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         cs = args[i];
         n = strlen(cs);
         for (k=0; k < n; k++)
         {
            switch(cs[k])
            {
            case 'm':
               *C |= 1<<CON_NOMVEC;
               break;
            case 'k':
               *C |= 1<<CON_NOKVEC;
               break;
            case 'c':
               *C |= 1<<CON_NOCOMPK;
               break;
            default:
               PrintUsage(args[0], -i, "UNKNOWN CONSTRAINT");
            }
         }
      case 'B':
         if (nbs)
         {
            free(mbs);
            free(nbs);
            free(kbs);
         }
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         N = atoi(args[i]);
         assert(N > 0);
         nbs = malloc(N*sizeof(int));
         assert(nbs);
         mbs = malloc(N*sizeof(int));
         assert(mbs);
         kbs = malloc(N*sizeof(int));
         assert(kbs);
         for (k=0; k < N; k++)
         {
           if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
           mbs[k] = atoi(args[i]);
           if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
           nbs[k] = atoi(args[i]);
           if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
           kbs[k] = atoi(args[i]);
         }
         break;
      case 'b':
         if (nbs)
         {
            free(mbs);
            free(nbs);
            free(kbs);
         }
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         N = atoi(args[i]);
         if (N < 0)
         {
            ALL = (N == -2) ? 32 : -1;  /* -2: don't force MB */
            N = 36;
         }
         assert(N > 0);
         nbs = malloc(N*sizeof(int));
         assert(nbs);
         mbs = malloc(N*sizeof(int));
         assert(nbs);
         kbs = malloc(N*sizeof(int));
         assert(nbs);
         if (ALL)  /* special case to just try most possible NBs */
         {
            mbs[ 0]=nbs[ 0]=kbs[ 0]=4;
            mbs[ 1]=nbs[ 1]=kbs[ 1]=6;
            mbs[ 2]=nbs[ 2]=kbs[ 2]=8;
            mbs[ 3]=nbs[ 3]=kbs[ 3]=12;
            mbs[ 4]=nbs[ 4]=kbs[ 4]=14;
            mbs[ 5]=nbs[ 5]=kbs[ 5]=16;
            mbs[ 6]=nbs[ 6]=kbs[ 6]=20;
            mbs[ 7]=nbs[ 7]=kbs[ 7]=22;
            mbs[ 8]=nbs[ 8]=kbs[ 8]=24;
            mbs[ 9]=nbs[ 9]=kbs[ 9]=26;
            mbs[10]=nbs[10]=kbs[10]=28;
            mbs[11]=nbs[11]=kbs[11]=32;
            mbs[12]=nbs[12]=kbs[12]=36;
            mbs[13]=nbs[13]=kbs[13]=40;
            mbs[14]=nbs[14]=kbs[14]=44;
            mbs[15]=nbs[15]=kbs[15]=48;
            mbs[16]=nbs[16]=kbs[16]=52;
            mbs[17]=nbs[17]=kbs[17]=56;
            mbs[18]=nbs[18]=kbs[18]=60;
            mbs[19]=nbs[19]=kbs[19]=64;
            mbs[20]=nbs[20]=kbs[20]=72;
            mbs[21]=nbs[21]=kbs[21]=80;
            mbs[22]=nbs[22]=kbs[22]=84;
            mbs[23]=nbs[23]=kbs[23]=88;
            mbs[24]=nbs[24]=kbs[24]=96;
            mbs[25]=nbs[25]=kbs[25]=104;
            mbs[26]=nbs[26]=kbs[26]=112;
            mbs[27]=nbs[27]=kbs[27]=120;
            mbs[28]=nbs[28]=kbs[28]=128;
            mbs[29]=nbs[29]=kbs[29]=132;
            mbs[30]=nbs[30]=kbs[30]=144;
            mbs[31]=nbs[31]=kbs[31]=156;
            mbs[32]=nbs[32]=kbs[32]=168;
            mbs[33]=nbs[33]=kbs[33]=192;
            mbs[34]=nbs[34]=kbs[34]=216;
            mbs[35]=nbs[35]=kbs[35]=240;
            for (k=0; k < ALL; k++)     /* don't force any particular */
               mbs[k] = 0;              /* MB during search */
         }
         else
         {
            for (k=0; k < N; k++)
            {
              if (++i >= nargs)
                  PrintUsage(args[0], i-1, NULL);
               mbs[k] = kbs[k] = nbs[k] = atoi(args[i]);
            }
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (B0 && mmb)
   {
      ATL_mmnode_t *mp;
      mp = BestBlocking_BFI(1, *PRE, mmb, B0, BN, KI, 0);
      KillMMNode(mmb);
      KillMMNode(mp);
      exit(0);
   }
   if (!nbs && !mmb)
      PrintUsage(args[0], -1, "Dimensional flag (-b or -B) or -S required!");
   *MBS = mbs;
   *NBS = nbs;
   *KBS = kbs;
   if (*PRE == 's' || *PRE == 'c')
   {
      VLEN[VTAVXZ] = 16;
      VLEN[VTAVX] = 8;
      VLEN[VTSSE] = 4;
      VLEN[VTGV] = 4;
      TSIZE = 4;
   }
   #ifdef ATL_AVXZ
      VECi = VTAVXZ;
   #elif defined(ATL_AVX)
      VECi = VTAVX;
   #elif defined(ATL_SSE1)
      if (*PRE == 's')
         VECi = VTSSE;
      #ifdef ATL_SSE2
      else
         VECi = VTSSE;
      #endif
   #elif (defined(ATL_AltiVec) && !defined(ATL_VSX)) || \
         (defined(ATL_NEON) && defined(ATL_NONIEEE) && ATL_NONIEEE != 0) || \
         (defined(ATL_3DNow) && defined(ATL_NONIEEE) && ATL_NONIEEE != 0)
      if (pre == 's')
         VECi = VTGV;
   #elif defined(ATL_VSX)
      VECi = VTGV;
   #endif
   *NN = N;
   return(mmb);
}

/*
 * Finds best-performing square cases from input list of good kernels, which
 * are first pruned of any losers (any kernel that is not faster than a
 * kernel with a smaller NB).
 */
ATL_mmnode_t *FindBestSquareKern(char pre, int verb, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mmp, *mmSQ=NULL, *prev=NULL;

   mmb = HarshPrune(mmb);  /* kill uncompetitive kernels */
   assert(mmb);

   mmSQ = FindBestSquareCases(pre, verb, mmb);
   assert(mmSQ);
   KillAllMMNodes(mmb);
   WriteMMFileWithPath(pre, "res", "uAMMFRC.sum", mmSQ);
   printf("DONE.\n");
   return(mmSQ);
}

int main(int nargs, char **args)
{
   int *nbs, *mbs, *kbs, verb, N, C=0, KCLEAN;
   float tol;
   char pre;
   ATL_mmnode_t *mmb, *mp, *mmB;
   double mfB;

   mmb = GetFlags(nargs, args, &pre, &verb, &N, &C, &mbs, &nbs, &kbs, &tol,
                  &KCLEAN);
/*
 * If -S is thrown, we take kernels from prior search, and just retune with
 * requirement that the kernels be square
 */
   if (mmb)
   {
      mmb = ApplyConstraints(mmb, C);
      assert(mmb);
      ApplyMoves2Flags(mmb, IMVS);
      mmb = FindBestSquareKern(pre, verb, mmb);
   }
   else
   {
      char upr=pre;
      if (pre == 'z')
         upr = 'd';
      else if (pre == 'c')
         upr = 's';
      mmb = GetWorkingUserCases(verb, upr);  /* get all working kernels */
      mmb = ApplyConstraints(mmb, C);        /* reject disallowed kerns */
      ApplyMoves2Flags(mmb, IMVS);           /* indicate A/B/C movememnt */
      mmb = FindBestKerns(pre, verb, mmb, N, mbs, nbs, kbs, C, tol);
   }
   if (KCLEAN)
      ComputeKClean(verb, pre, mmb);
   mmb = ATL_SortMMNodesByMflop(0, mmb);
   mmB = ReadMMFileWithPath(pre, "res", "eAMMRES.sum");
   mmB = ATL_SortMMNodesByMflop(0, mmB);
   printf("\n");
   if (mmb)
      printf("YOUR BEST CASE: %d-%s, BLK=%d,%d,%d  %.2f\n",
             mmb->ID, mmb->rout, mmb->mbB, mmb->nbB, mmb->kbB, mmb->mflop[0]);
   if (mmB)
   {
      printf("ATLAS BEST CASE: %d-%s, BLK=%d,%d,%d  %.2f\n", mmB->ID,
             mmB->rout, mmB->mbB, mmB->nbB, mmB->kbB, mmB->mflop[0]);
      if (mmb)
         printf("   RATIO=%0.4f\n", mmb->mflop[0] / mmB->mflop[0]);
   }

   KillAllMMNodes(mmb);
   KillAllMMNodes(mmB);
   free(mbs);
   free(nbs);
   free(kbs);
   return(0);
}
