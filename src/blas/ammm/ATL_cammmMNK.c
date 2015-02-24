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
 * This routine called when N < M and K is large
 */
int Mjoin(PATL,ammmMNK)
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
   size_t m, nsmblks, nmblks, nsnblks, nnblks, i, incAm0, incAm, incAw0;
   size_t nkb, incAk, incAk0, mulAm, incBk, incBk0, mulBn;
   int mu, nu, ku, MB, NB, KB, KB0, kb0, NMU, NNU, A_1TRIP;
   void *vp;
   TYPE *rC, *iC, *iB, *pB0, *iA, *pA0;
   ammkern_t ammK0, ammK0_bn, ammK0_b1, amm_b1, amm_bn;
   const int B_BYCOLS = (TB == AtlasNoTrans || TB == AtlasConj);
   const int A_BYROWS = (TA == AtlasNoTrans || TA == AtlasConj);
   const TYPE ONE[2] = {ATL_rone, ATL_rzero};
   const TYPE *alpA=ONE, *alpB=ONE, *alpC=ONE;
   amminfo_t mminfo;

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
   NB = ATL_ComputeB(N, nu, MY_MAXNB, &nsnblks, &nnblks);
   NNU = NB / nu;
   KB = mminfo.kb;
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
            ammK0 = mminfo.amm_b0;
         else
            ammK0 = (mminfo.kb==KB0) ? mminfo.amm_b0 : mminfo.amm_k1_b0;
      }
      else
   #endif
   {
      ammK0 = (kb0 == KB) ? mminfo.amm_b0 : mminfo.amm_k1_b0;
      if (ATL_AMMFLG_KRUNTIME(mminfo.flag) && kb0 == (kb0/ku)*ku &&
          kb0 > mminfo.kbmin)
         ammK0 = mminfo.amm_b0;
   }
   if (ammK0 == mminfo.amm_b0)
   {
      amm_b1 = ammK0_b1 = mminfo.amm_b1;
      amm_bn = ammK0_bn = mminfo.amm_bn;
   }
   else
   {
      ammK0_b1 = mminfo.amm_k1_b1;
      ammK0_bn = mminfo.amm_k1_bn;
      amm_b1 = mminfo.amm_b1;
      amm_bn = mminfo.amm_bn;
   }
   a2blk = mminfo.a2blk;
   b2blk = mminfo.b2blk;
   blk2c = mminfo.Cblk2cm;
   i = nkb*KB+KB0;
/*
 * Determine worspace requirements and allocate
 */
   {
      size_t tsz;
      const size_t szA=MB*i;
      const size_t szB=i*(nsnblks*(NB-nu)+(nnblks-nsnblks)*NB);
      const int szC = MB*NB;

      tsz = ATL_MulBySize(szA + szB + szC + mu*nu*ku) + 3*ATL_Cachelen;
      if (tsz > ATL_MaxMalloc)
         return(2);
      vp = malloc(tsz);
      if (!vp)
         return(1);
      iC = ATL_AlignPtr(vp);
      rC = iC + szC;
      iA = rC + szC;
      pA0 = iA = ATL_AlignPtr(iA);
      iB = iA + szA + szA;
      pB0 = iB = ATL_AlignPtr(iB);
   }
   if (A_BYROWS)
   {
      incAk = KB*(lda SHIFT);
      incAk0 = kb0*(lda SHIFT);
      mulAm = 1 SHIFT;
   }
   else
   {
      incAk = KB SHIFT;
      incAk0 = kb0 SHIFT;
      mulAm = lda SHIFT;
   }
   if (B_BYCOLS)
   {
      incBk = KB SHIFT;
      incBk0 = kb0 SHIFT;
      mulBn = ldb SHIFT;
   }
   else
   {
      incBk = (KB SHIFT)*ldb;
      incBk0 = kb0*(ldb SHIFT);
      mulBn = 1 SHIFT;
   }

   for (m=M, i=0; i < nmblks; i++)
   {
      size_t j, n;
      int mb, mm, nmu, incAw, incAw0;
      TYPE *c=C;
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
      mm = Mmin(m, mb);  /* number of A/C rows left */
      m -= mm;
      incAw = mb*KB;
      incAw0 = mb*KB0;
      for (n=N, j=0; j < nnblks; j++)
      {
         size_t k;
         int nb, nn, nnu, incBw, incBw0;
         const TYPE *b=B, *a=A;
         TYPE *pAn, *pBn, *rA, *rB;

         if (j < nsnblks)
         {
            nb = NB-nu;
            nnu = NNU-1;
         }
         else
         {
            nb = NB;
            nnu = NNU;
         }
         incBw = KB*nb;
         nn = Mmin(n, nb);  /* number of B/C cols left */
         n -= nn;
         incBw0 = KB0*nb;

         rA = iA + incAw0;
         pAn = rA + incAw0;
         pAn = (nkb) ? pAn : pA0;
         rB = iB + incBw0;
         pBn = rB + incBw0;
         if (!j)
         {
            a2blk(kb0, mm, alpA, a, lda, rA, iA);
            a += incAk0;
         }
         if (!i)
         {
             b2blk(kb0, nn, alpB, b, ldb, rB, iB);
             b += incBk0;
         }
         ammK0(nmu, nnu, KB0, iA, iB, rC, rA, iB, iC);
         ammK0(nmu, nnu, KB0, rA, iB, iC, rA, rB, rC);
         ammK0_bn(nmu, nnu, KB0, rA, rB, rC, iA, rB, iC);
         ammK0_b1(nmu, nnu, KB0, iA, rB, iC, pAn, pBn, rC);
         iA = pAn;
         iB = pBn;
         for (k=0; k < nkb; k++)
         {
            rA = iA + incAw;
            rB = iB + incBw;
            if (!j)
            {
               a2blk(KB, mm, alpA, a, lda, rA, iA);
               a += incAk;
            }
            if (!i)
            {
                b2blk(KB, nn, alpB, b, ldb, rB, iB);
                b += incBk;
            }
            pAn = rA + incAw;
            pAn = (k != nkb-1) ? pAn : pA0;
            pBn = rB + incBw;
            pBn = (k != nkb-1 || j != nnblks-1) ? pBn : pB0;
            amm_bn(nmu, nnu, KB, iA, iB, rC, rA, iB, iC);
            amm_b1(nmu, nnu, KB, rA, iB, iC, rA, rB, rC);
            amm_bn(nmu, nnu, KB, rA, rB, rC, iA, rB, iC);
            amm_b1(nmu, nnu, KB, iA, rB, iC, pAn, pBn, rC);
            iA = pAn;
            iB = pBn;
         }
         blk2c(mm, nn, alpC, rC, iC, beta, c, ldc);
         c += nn*(ldc SHIFT);
         B += nn*mulBn;
      }
      iB = pB0;
      A += mm*mulAm;
      C += mm SHIFT;
   }

   free(vp);
   return(0);
}
