#include "atlas_misc.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_sum.h))
/*
 * This routine called in degenerate case where all dims less than max block,
 * so we can do entire operation with one kernel call
 */
int Mjoin(PATL,ammm_1b)
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
   int i;
   int nmu, nnu, nku, bM, bN, bK, szB;
   #if ATL_AMM_MAXKMAJ > 1
      int KK;
   #else
      #define KK K
   #endif
   int mu, nu, ku, appAl;
   void *vp;
   TYPE *pA, *pB, *pC, *p;
   TYPE alpA=ATL_rone, alpB=ATL_rone, alpC=ATL_rone;
   ammkern_t amm;
   ablk2cmat_t blk2c;
   cm2am_t a2blk, b2blk;
   amminfo_t mminfo;
   appAl = Mjoin(PATL,GetAmmmInfo)(&mminfo, TA, TB, M, N, K, alpha, beta);
   if (!appAl)
      alpA = alpha;
   else if (appAl == 1)
      alpB = alpha;
   else
      alpC = alpha;
/*
 * These kernels all take runtime M/N, and do well with near-square, so
 * blindly use this kernel with nM = CEIL(M/mu)*mu, nN = CEIL(N/nu)*nu
 */
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   nmu = (M+mu-1)/mu;
   nnu = (N+nu-1)/nu;
   nku = (K+ku-1)/ku;
   bM = nmu * mu;
   bN = nnu * nu;
   bK = nku * ku;
   a2blk = mminfo.a2blk;
   b2blk = mminfo.b2blk;
   blk2c = mminfo.Cblk2cm;
   #if ATL_AMM_MAXKMAJ > 1
      KK = K;
      if (ATL_AMMFLG_KMAJOR(mminfo.flag))
      {
         KK = bK;
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag))
            amm = mminfo.amm_b0;
         else
            amm = (mminfo.kb==KK) ? mminfo.amm_b0 : mminfo.amm_k1_b0;
      }
      else
   #endif
         amm = (bK == K) ? mminfo.amm_b0 : mminfo.amm_k1_b0;
/*
 * Force rank-K code to handle this case if we would have to use K-cleanup
 * code with unknown performance
 */
   if (amm != mminfo.amm_b0 && K > 2 && K <= ATL_MAXK_RKK)
      return(-1);
   szB = KK*bN;
   vp = malloc((bM*bN + bM*KK + szB + mu*nu*ku)*sizeof(TYPE) + 3*ATL_Cachelen);
   if (!vp)
      return(1);
   pB = ATL_AlignPtr(vp);
   pA = pB + szB;
   pA = ATL_AlignPtr(pA);
   pC = pA + bM*KK;
   pC = ATL_AlignPtr(pC);
/*
 * Copy A & B into workspace, and pad K its if necessary
 */
   a2blk(K, M, alpA, A, lda, pA);
   b2blk(K, N, alpB, B, ldb, pB);
   amm(nmu, nnu, KK, pA, pB, pC, pA, pB, pC);
   blk2c(M, N, alpC, pC, beta, C, ldc);

   free(vp);
   return(0);
}
