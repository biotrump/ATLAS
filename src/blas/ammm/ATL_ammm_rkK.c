#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_level1.h"
#include "atlas_lvl2.h"
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_sum.h))

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
 * This routine called when 2 < K <= MAXK
 */
int Mjoin(PATL,ammm_rkK)
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
   size_t n, nsmblks, nmblks, nsnblks, nnblks, j, incAm0, incAm, incAw0, incAw;
   int mu, nu, ku, MB, NB, mb, nb, NMU, NNU, A_1TRIP;
   void *vp;
   TYPE *pC, *pB, *pA, *pA0;
   ammkern_t amm;
   const int B_BYCOLS = (TB == AtlasNoTrans);
   const int A_BYROWS = (TA == AtlasNoTrans);
   TYPE alpA=ATL_rone, alpB=ATL_rone;
   #if ATL_MAXKMAJ_RKK > 1
      int KK=K;
   #else
      #define KK K
   #endif
   amminfo_t mminfo;

   mu = Mjoin(PATL,GetRankKInfo)(&mminfo, TA, TB, M, N, K, alpha, beta);
   if (!mu)
      alpA = alpha;
   else
      alpB = alpha;
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   MB = ATL_ComputeB(M, mu, ATL_MAXM_RKK, &nsmblks, &nmblks);
   NMU = MB / mu;
   NB = ATL_ComputeB(N, nu, ATL_MAXN_RKK, &nsnblks, &nnblks);
   NNU = NB / nu;
   A_1TRIP = (nnblks < 2);
   #if ATL_MAXKMAJ_RKK > 1
      if (ATL_AMMFLG_KMAJOR(mminfo.flag))
         KK = ((K+ku-1)/ku)*ku;
   #endif

   {
      size_t tsz, szA;
      int szB = KK*NB, szC = MB*NB;

      if (A_1TRIP)
         szA = MB*KK;
      else
         szA = (nsmblks*(MB-mu)+(nmblks-nsmblks)*MB)*KK;
      tsz = ATL_MulBySize(szA + szB + szC + mu*nu*ku) + 3*ATL_Cachelen;
      if (tsz > ATL_MaxMalloc)
         return(2);
      vp = malloc(tsz);
      if (!vp)
         return(1);
      pC = ATL_AlignPtr(vp);
      pB = pC + szC;
      pB = ATL_AlignPtr(pB);
      pA = pB + szB;
      pA0 = pA = ATL_AlignPtr(pA);
   }
   amm = mminfo.amm_b0;
   a2blk = mminfo.a2blk;
   b2blk = mminfo.b2blk;
   blk2c = mminfo.Cblk2cm;

   if (A_BYROWS)
   {
      incAm0 = (MB-mu);
      incAm = MB;
   }
   else
   {
      incAm0 = (MB-mu)*lda;
      incAm = MB*lda;
   }
   if (A_1TRIP)
      incAw0 = incAw = 0;
   else
   {
      incAw0 = (MB-mu)*KK;
      incAw = MB*KK;
   }

   C += N*ldc;
   B += (B_BYCOLS) ? N*ldb : N;
   n = N;
   for (n=N, j=0; j < nnblks; j++)
   {
      int nmu, mb, nnu, nb, nn;
      size_t incBn, i, m;
      TYPE *c;
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
      nn = Mmin(n, nb);
      B -= (B_BYCOLS) ? nn*ldb : nn;
      C -= nn*ldc;
      c=C;

      b2blk(K, nn, alpB, B, ldb, pB);

      mb = MB-mu;
      nmu = NMU-1;
      pA = pA0;
      for (m=M,i=0; i < nsmblks; i++, m -= mb, c += mb)
      {
         const int mm = Mmin(m, mb);
         TYPE *pAn = pA + incAw0;
         if (!j)
         {
            a2blk(K, mm, alpA, A, lda, pA);
            A += incAm0;
         }
         amm(nmu, nnu, KK, pA, pB, pC, pAn, pB, pC);
         blk2c(mm, nn, ATL_rone, pC, beta, c, ldc);
         pA = pAn;
      }

      for (mb=MB, nmu=NMU; i < nmblks; i++, m -= mb, c += mb)
      {
         const int mm = Mmin(m, mb);
         TYPE *pAn = pA + incAw;
         if (!j)
         {
            a2blk(K, mm, alpA, A, lda, pA);
            A += incAm;
         }
         amm(nmu, nnu, KK, pA, pB, pC, pAn, pB, pC);
         blk2c(mm, nn, ATL_rone, pC, beta, c, ldc);
         pA = pAn;
      }

      n -= nb;
   }
   free(vp);
   return(0);
}
static void *FixVector(enum ATLAS_TRANS TX, ATL_CSZT N, const SCALAR alpha,
                       const TYPE *X, ATL_CSZT incX)
{
   void *vx;
   TYPE *x;
   vx = malloc(ATL_MulBySize(N)+ATL_Cachelen);
   ATL_assert(vx);
   x = ATL_AlignPtr(vx);
   if (SCALAR_IS_ONE(alpha))
      Mjoin(PATL,copy)(N, X, incX, x, 1);
   else
      Mjoin(PATL,cpsc)(N, alpha, X, incX, x, 1);
   return(vx);
}

/*
 * This entry makes rkK safe for L3kernel aliased calls.  It handles
 * only the aliasing required by the L3kernels, namely square blocks
 * less than ATLAS's largest blocking factor for the square dimensions,
 * with one of A/B aliased with C, and aliased by having
 * either A == C or B == C (i.e., not a partial overlap).  When A==C,
 * M=K < ATL_MAXK_RKK; when B==C, N=K < ATL_MAXK_RKK.
 */
int Mjoin(PATL,ammm_aliased_rkK)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alp,
   const TYPE *A,
   ATL_CSZT lda,
   const TYPE *B,
   ATL_CSZT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   void *vp=NULL;

   if (K == 0 || SCALAR_IS_ZERO(alp))
   {
      if (SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gezero)(M, N, C, ldc);
      else if (!SCALAR_IS_ZERO(beta))
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);
      return(0);
   }
   if (K == 1)
   {
      TYPE alpha = alp;
      void *xp=NULL, *yp=NULL;
      TYPE *a = (TYPE*)A, *b = (TYPE*)B;
      size_t ldA=lda, ldB=ldb;


      if (A == C)
      {
         xp = FixVector(TA, M, alpha, A,
                        (TA == AtlasTrans || TA == AtlasConjTrans) ? lda:1);
         a = ATL_AlignPtr(xp);
         alpha = ATL_rone;
         if (TA == AtlasConjTrans || TA == AtlasTrans)
         {
            TA = AtlasTrans;
            ldA = 1;
         }
         else if (TA == AtlasConj)
            TA = AtlasNoTrans;
      }
      if (B == C)
      {
         yp = FixVector(TB, N, alpha, B,
                        (TB == AtlasTrans || TB == AtlasConjTrans) ? 1:ldb);
         b = ATL_AlignPtr(yp);
         alpha = ATL_rone;
         if (TB == AtlasConjTrans)
            TB = AtlasTrans;
         else if (TB == AtlasConj || TB == AtlasNoTrans)
         {
            TA = AtlasNoTrans;
            ldB = 1;
         }
      }
      ATL_assert(!Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, a, ldA, b, ldB,
                                   beta, C, ldc));
      if (xp)
         free(xp);
      if (yp)
         free(yp);
      return(0);
   }

   if (K == 2)
   {
      TYPE alpha = alp;
/*
 *    If BETA != 1, ammm_rk2 will copy all inputs and thus aliasing safe
 */
      if (!SCALAR_IS_ONE(beta))
         Mjoin(PATL,ammm_rk2)(TA, TB, M, N, alpha, A, lda, B, ldb, beta,
                              C, ldc);
/*
 *    For beta = 1, copy aliased input array(s) and then call GER2
 */
      else
      {
         void *wp=NULL, *xp=NULL, *yp=NULL, *zp=NULL;
         TYPE *w, *x, *y, *z;
         ATL_SZT incX, incY;
         if (A == C)
         {
            ATL_CSZT incx = (TA == AtlasTrans) ? lda:1;
            wp = FixVector(TA, M, alpha, A, incx);
            xp = FixVector(TA, M, alpha, A+((incx==1 ? lda:1)SHIFT), incx);
            alpha = ATL_rone;
            w = ATL_AlignPtr(wp);
            incX = 1;
         }
         else  /* don't need to copy A */
         {
            w = (TYPE*)A;
            if (TA == AtlasNoTrans)
            {
               x = (TYPE*)(A + (lda SHIFT));
               incX = 1;
            }
            else /* if (TA == AtlasTrans) */
            {
               x = (TYPE*)(A + (1 SHIFT));
               incX = lda;
            }
         }
         if (B == C)
         {
            ATL_CSZT incy = (TB == AtlasTrans || TB == AtlasConjTrans) ? 1:ldb;
            yp = FixVector(TB, N, alpha, A, incy);
            zp = FixVector(TB, N, alpha, A+(((incy==1)?ldb:1)SHIFT), incy);
            y = ATL_AlignPtr(yp);
            z = ATL_AlignPtr(zp);
            incY = 1;
            alpha = ATL_rone;

         }
         else  /* no need to copy B */
         {
            y = (TYPE*)B;
            if (TB == AtlasNoTrans)
            {
               incY = ldb;
               z = (TYPE*)(B + (1 SHIFT));
            }
            else
            {
               incY = 1;
               z = (TYPE*)(B + (ldb SHIFT));
            }
         }
         Mjoin(PATL,ger2)(M, N, alpha, w, incX, y, incY, ATL_rone,
                          x, incX, z, incY, C, ldc);
         if (wp)
            free(wp);
         if (xp)
            free(xp);
         if (yp)
            free(yp);
         if (zp)
            free(zp);
         return(0);
      }
      return(0);
   }
/*
 * For K > 3, ATL_ammm_rkK is safe for these precise aliasing conditions
 */
   ATL_assert(!Mjoin(PATL,ammm_rkK)(TA, TB, M, N, K, alp, A, lda, B, ldb,
                                    beta, C, ldc));
   return(0);
}
