#include "atlas_misc.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_sum.h))
/*
 * This routine handles M <= MAXM && N <= MAXN && very long K, or
 * the inner-product GEMM form.  It appears in the GEMM-based SYRK, which
 * is important for Cholesky.  It is typically the worst-case for ATLAS,
 * since the copy of A and B are of the same order as the computation.
 * It is a minimal workspace routine.
 */
int Mjoin(PATL,ammm_IP)
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
   ammkern_t ammK0, amm;
   amminfo_t mminfo;
   TYPE alpA=ATL_rone, alpB=ATL_rone, alpC=ATL_rone;
   TYPE *pA, *pB, *pC;
   int mu, nu, ku, nmu, nnu, MM, NN, KB, KB0, kb0;
   size_t incA, incB, nkb, k;
   void *vp;

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
   nmu = (M+mu-1)/mu;
   nnu = (N+nu-1)/nu;
   KB = mminfo.kb;
   MM = nmu * mu;
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
   amm = mminfo.amm_b1;
   a2blk = mminfo.a2blk;
   b2blk = mminfo.b2blk;
   blk2c = mminfo.Cblk2cm;

/*
 * Do memory allocation, setup pointers
 */
   {
      const int szA=MM*KB, szB=KB*NN, szC=MM*NN;
      vp = malloc(ATL_MulBySize(szA + szB + szC + mu*nu*ku) + 3*ATL_Cachelen);
      ATL_assert(vp);
      pA = ATL_AlignPtr(vp);
      pB = pA + szA;
      pB = ATL_AlignPtr(pB);
      pC = pB + szB;
      pC = ATL_AlignPtr(pC);
   }
   incA = (TA == AtlasNoTrans) ? lda*KB : KB;
   incB = (TB == AtlasNoTrans) ? KB : KB*ldb;
/*
 * Do first (possibly partial) K-block
 */
   a2blk(kb0, M, alpA, A, lda, pA);
   b2blk(kb0, N, alpB, B, ldb, pB);
   ammK0(nmu, nnu, KB0, pA, pB, pC, pA, pB, pC);
   A += (TA == AtlasNoTrans) ? lda*kb0 : kb0;
   B += (TB == AtlasNoTrans) ? kb0 : kb0*ldb;
/*
 * Loop over all full-sized blocks
 */
   for (k=0; k < nkb; k++)
   {
      a2blk(KB, M, alpA, A, lda, pA);
      b2blk(KB, N, alpB, B, ldb, pB);
      amm(nmu, nnu, KB, pA, pB, pC, pA, pB, pC);
      A += incA;
      B += incB;
   }
   blk2c(M, N, alpC, pC, beta, C, ldc);

   free(vp);
   return(0);
}
