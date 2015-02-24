#include "atlas_misc.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_sum.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_sum.h))
#ifdef ATL_CAMM_MAXMB
   #define MY_MAXMB ATL_CAMM_MAXMB
   #define MY_MAXNB ATL_CAMM_MAXNB
#else
   #define MY_MAXMB ATL_MAXM_RKK
   #define MY_MAXNB ATL_MAXM_RKK
#endif

#ifdef __GNUC__
static inline int ATL_ComputeB   /* RETURNS: selected blocking */
#else
static int ATL_ComputeB           /* RETURNS: selected blocking */
#endif
(
   size_t N,   /* problem dimension */
   int nu,     /* unrolling by kernel on this dim */
   int nb,     /* IN: large-case blocking */
   size_t *NS, /* OUT: # of blks of size NB-nu to perform */
   size_t *NT  /* OUT: # of blks to perform */
)
{
   size_t ns, nt, nblks, NN;
/*
 * If the entire problem is less than or equal to the unrolling, choose a block
 * of the ceiling of the unrolling and only do one
 */
   NN=((N+nu-1)/nu)*nu;  /* ceiling of number of unrollings in N */
   if (NN <= nu)
   {
      *NS = 0;
      *NT = 1;
      return(NN);
   }
/*
 * If suggested block size is smaller or same as unrolling, then the blocking
 * size is the unrolling, and we don't have an NB-nu sized-blocks, since that
 * would be zero sized
 */
   if (nb <= nu)
   {
      *NS = 0;
      *NT = NN/nu;
      return(nu);
   }

   nb = (nb/nu)*nu;      /* floor of number of unrollings in a block*/
/*
 * If 1 block is within NU of covering the entire dim, just make the
 * block size the entire dim
 */
   if (nb+nu >= NN)
   {
      *NS = 0;
      *NT = 1;
      return(NN);
   }
/*
 * Otherwise, compute how many blocks we need  of each type
 */
   while(1)
   {
      nblks = (N+nb-1)/nb;
      ns = (nblks*nb - NN)/nu;
      if (ns < nblks)
         break;
      nb -= nu;
   }

   *NS = ns;
   *NT = nblks;
   return(nb);
}
/*
 * This routine handles N <= MAXN, K & M large (left-looking shape)
 */
int Mjoin(PATL,ammm_tN)
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
   ablk2cmat_t blk2c;
   cm2am_t a2blk, b2blk;
   ammkern_t ammK0_b0, ammK0_b1, ammK0_bn, amm_b1, amm_bn;
   amminfo_t mminfo;
   const TYPE ONE[2] = {ATL_rone, ATL_rzero};
   const TYPE *alpA=ONE, *alpB=ONE, *alpC=ONE;
   TYPE *rA, *iA, *pB0, *rC, *iC;
   int mu, nu, ku, nnu, NN, MB, NMU, KB, KB0, kb0, incBw, incBw0;
   size_t incAk0, incAk, mulAm, incBk0, incBk, nkb, k, nmblks, nsmblks, i, m;
   void *vp;

   ATL_assert(N <= ATL_AMM_MAXNB);
   mu = Mjoin(PATL,GetAmmmInfo)(&mminfo, TA, TB, M, N, K, alpha, beta);
   if (!mu)
      alpA = alpha;
   else if (mu == 1)
      alpB = alpha;
   else
      alpC = alpha;

   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   MB = ATL_ComputeB(M, mu, MY_MAXMB, &nsmblks, &nmblks);
   NMU = MB / mu;
   nnu = (N+nu-1)/nu;
   KB = mminfo.kb;
   NN = nnu * nu;
   nkb = K/KB;
/*
 * kb0: K remainder, KB0 is CEIL(kb0/ku)*ku for k-vector kerns, and
 * same as kb0 for M-vector kerns
 */
   KB0 = kb0 = K - nkb*KB;
   if (!kb0)
   {
      KB0 = kb0 = KB;
      nkb--;
   }
   #if ATL_AMM_MAXKMAJ > 1
      if (ATL_AMMFLG_KMAJOR(mminfo.flag))
      {
         KB0 = ((kb0+ku-1)/ku)*ku;
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag))
            ammK0_b0 = mminfo.amm_b0;
         else
            ammK0_b0 = (mminfo.kb==KB0) ? mminfo.amm_b0 : mminfo.amm_k1_b0;
      }
      else
   #endif
   {
      ammK0_b0 = (kb0 == KB) ? mminfo.amm_b0 : mminfo.amm_k1_b0;
      if (ATL_AMMFLG_KRUNTIME(mminfo.flag) && (kb0/ku)*ku == kb0 &&
          kb0 > mminfo.kbmin)
         ammK0_b0 = mminfo.amm_b0;
   }
   amm_b1 = mminfo.amm_b1;
   amm_bn = mminfo.amm_bn;
   if (ammK0_b0 == mminfo.amm_b0)
   {
      ammK0_b1 = mminfo.amm_b1;
      ammK0_bn = mminfo.amm_bn;
   }
   else
   {
      ammK0_b1 = mminfo.amm_k1_b1;
      ammK0_bn = mminfo.amm_k1_bn;
   }
   a2blk = mminfo.a2blk;
   b2blk = mminfo.b2blk;
   blk2c = mminfo.Cblk2cm;

/*
 * Do memory allocation, setup pointers
 */
   {
      const int szA = MB*KB, szC=MB*NN;
      size_t szB = (nkb*KB+KB0)*NN;
      const size_t tsz = ATL_MulBySize(szA+szB+szC+mu*nu*ku) + 3*ATL_Cachelen;
      if (tsz > ATL_MaxMalloc)
         return(2);
      vp = malloc(tsz);
      if (!vp)
         return(1);
      iA = ATL_AlignPtr(vp);
      rA = iA + szA;
      iC = rA + szA;
      iC = ATL_AlignPtr(iC);
      rC = iC + szC;
      pB0 = rC + szC;
      pB0 = ATL_AlignPtr(pB0);
   }
   if (TA == AtlasNoTrans)
   {
      incAk = (lda*KB)SHIFT;
      incAk0 = (lda*kb0)SHIFT;
      mulAm = 2;
   }
   else
   {
      incAk = KB SHIFT;
      incAk0 = kb0 SHIFT;
      mulAm = lda SHIFT;
   }
   if (TB == AtlasNoTrans)
   {
      incBk = KB SHIFT;
      incBk0 = kb0 SHIFT;
   }
   else
   {
      incBk = (KB*ldb)SHIFT;
      incBk0 = (kb0*ldb)SHIFT;
   }
   incBw0 = KB0*NN;
   incBw = KB*NN;
   for (m=M,i=0; i < nmblks; i++)
   {
      const TYPE *An;
      TYPE *iB=pB0, *rB=pB0+incBw0, *pBn=rB+incBw0;
      int mb, nmu, mm;

      if (i < nsmblks)
      {
         mb = MB-mu;
         nmu = NMU-1;
      }
      else
      {
         mb = MB;
         nmu = NMU;
      }
      mm = Mmin(m, mb);
      m -= mm;
      An = A + mm*mulAm;
/*
 *    Do first (possibly partial) K-block
 */
      a2blk(kb0, mm, alpA, A, lda, rA, iA);
      A += incAk0;
      if (!i)
      {
         b2blk(kb0, N, alpB, B, ldb, rB, iB);
         B += incBk0;
      }
      rB = iB + incBw0;
      pBn = rB + incBw0;
      ammK0_b0(nmu, nnu, KB0, iA, iB, rC, rA, iB, iC);
      ammK0_b0(nmu, nnu, KB0, rA, iB, iC, rA, rB, rC);
      ammK0_bn(nmu, nnu, KB0, rA, rB, rC, iA, rB, iC);
      ammK0_b1(nmu, nnu, KB0, iA, rB, iC, rA, pBn, iC);
      iB = pBn;
/*
 *    Loop over all full-sized blocks
 */
      for (k=0; k < nkb; k++)
      {
         a2blk(KB, mm, alpA, A, lda, rA, iA);
         A += incAk;
         rB = iB + incBw;
         if (!i)
         {
            b2blk(KB, N, alpB, B, ldb, rB, iB);
            B += incBk;
         }
         pBn = (k < nkb-1) ? rB+incBw : pB0;
         amm_bn(nmu, nnu, KB, iA, iB, rC, rA, iB, iC);
         amm_b1(nmu, nnu, KB, rA, iB, iC, rA, rB, rC);
         amm_bn(nmu, nnu, KB, rA, rB, rC, iA, rB, iC);
         amm_b1(nmu, nnu, KB, iA, rB, iC, rA, pBn, iC);
         iB = pBn;
      }
      blk2c(mm, N, alpC, rC, iC, beta, C, ldc);
      C += mm+mm;
      A = An;
   }

   free(vp);
   return(0);
}
