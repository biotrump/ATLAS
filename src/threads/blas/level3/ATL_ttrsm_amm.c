#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"

static void cpyAblk(ATL_ttrsm_amm_t *pd)
{
   const int mb=pd->mb;
   int actr, ablk, kb;
   size_t i, j;
   const TYPE *a=pd->A;
   TYPE *wA;

   actr = ATL_DecAtomicCount(pd->AblkCtr);
   if (!actr)
      return;
   ablk = pd->nablks - actr;
   wA = pd->wA + ablk*pd->blkszA;
   Mtblk2coord(pd->nmblks, ablk, i, j);
   i = pd->mb0 + (i-1)*mb;
   if (j)
   {
      j = pd->mb0 + (j-1)*mb;
      kb = mb;
   }
   else
      kb = pd->mb0;
   if (pd->uplo == AtlasLower)
      a += j*pd->lda + i;
   else
      a += i*pd->lda + j;
   pd->a2blk(kb, mb, ATL_rnone, a, pd->lda, wA);
   ATL_mutex_lock(pd->Acpymut);
   ATL_SetBitBV(pd->AcpyBV, ablk);
   if (ATL_FindFirstUnsetBitBV(pd->AcpyBV, 0) == -1)
      pd->AcpyDone = 1;
   ATL_mutex_unlock(pd->Acpymut);
}

#define trsmK Mjoin(PATL,trsm)
void Mjoin(PATL,DoTrsm_amm)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_ttrsm_amm_t *pd = pp->PD;
   TYPE *wrk;
   int rhs;
   void *rhsCtr = pd->rhsCtr;

   wrk = pd->w + vrank*pd->wsL;
   if ((rhs = ATL_DecGlobalAtomicCount(rhsCtr, vrank)))
   {
      TYPE *wB=wrk, *X = pd->X;
      TYPE alpha = pd->alpha;
      const size_t ldx=pd->ldx, lda=pd->lda;
      const int mb=pd->mb, mb0=pd->mb0, MB0=pd->MB0, nmu=pd->nmu;
      const int nmblks=pd->nmblks, nnblks=pd->nnblks;
      const int blkszC = pd->blkszB, blkszA = pd->blkszA;
      const enum ATLAS_DIAG diag = pd->diag;
      const enum ATLAS_DIAG uplo = pd->uplo;
      const enum ATLAS_TRANS TA=pd->TA;
      cm2am_t b2blk = pd->b2blk;
      ablk2cmat_t blk2c=pd->blk2c;
      ammkern_t amm_b0=pd->amm_b0, amm = pd->amm_b1;

      do
      {
         const int rblk = nnblks - rhs;
         const int nb = (rblk == nnblks-1) ? pd->nbf : pd->nb;
         const int nnu = (rblk == nnblks-1) ? pd->nnuf : pd->nnu;
         TYPE *x = X + ldx*rblk*pd->nb;
         TYPE *wA=pd->wA, *wC=wB+blkszC;
         const TYPE *A = pd->A;
/*
 *       Solve top block using first diag blk of A (no need to copy)
 */
         trsmK(AtlasLeft, uplo, TA, diag, mb0, nb, alpha, A, lda, x, ldx);
/*
 *       If there are more RHS blocks, subtract solved X from them
 */
         if (nmblks > 1)
         {
            int i;
            TYPE *wc=wC, *a = wA;
            int ba=0;
/*
 *          Copy solved part of X, and subtract it from C, which will later be
 *          used to update the B values before solving them to X
 */
            b2blk(mb0, nb, ATL_rone, x, ldx, wB);
            for (i=1; i < nmblks; i++, ba++)
            {
               TYPE *cn = (i != nmblks-1) ? wc + blkszC:wC, *an = a + blkszA;
/*
 *             If the block of A we need hasn't been copied, copy A until it is.
 */
               if (!pd->AcpyDone)
               {
                  while(!ATL_IsBitSetBV(pd->AcpyBV, ba))
                     cpyAblk(pd);
               }
               amm_b0(nmu, nnu, MB0, a, wB, wc, an, wB, cn);
               wc = cn;
               a = an;
            }
/*
 *          For each remaining RHS block, we solve it, then subtract its
 *          part of the equation from all the blocks below it, until alg
 *          finishes.
 */
            A += mb0*(lda+1);
            x += mb0;
            for (i=1; i < nmblks; i++)
            {
/*
 *             Apply alpha to X, and then subtract of solved equations, before
 *             solving using the current mbxmb diagonal block of A to form X
 */
               blk2c(mb, nb, ATL_rone, wC, alpha, x, pd->ldx);
               wC += blkszC;
               trsmK(AtlasLeft, uplo, TA, diag, mb, nb, ATL_rone, A, lda,
                     x, ldx);
/*
 *             Subtract solved equations from remaining RHS blocks
 */
               if (i != nmblks-1)
               {
                  int k;

                  b2blk(mb, nb, ATL_rone, x, ldx, wB);
                  wc = wC;
                  for (k=i+1; k < nmblks; k++, ba++)
                  {
                     TYPE *cn = (k != nmblks-1) ? wc+blkszC : wC;
                     TYPE *an = a + blkszA;
                     if (!pd->AcpyDone)
                     {
                        while(!ATL_IsBitSetBV(pd->AcpyBV, ba))
                           cpyAblk(pd);
                     }
                     amm(nmu, nnu, mb, a, wB, wc, an, wB, cn);
                     wc = cn;
                     a = an;
                  }
               }
               A += mb*(lda+1);
               x += mb;
            }
         }
      }
      while ((rhs = ATL_DecGlobalAtomicCount(rhsCtr, vrank)));
   }
/*
 * See if any A remains to be copying before finishing; only if I never got
 * any work is it possible to reach here before A is completely copied.
 */
   else if (!pd->AcpyDone)
      while (ATL_GetAtomicCount(pd->AblkCtr))
         cpyAblk(pd);
}

/*
 * This routine handles case where A is on left of B, A X = B
 */
static int ttrsm_ammL
   (const enum ATLAS_UPLO uplo, const enum ATLAS_TRANS TA,
    const enum ATLAS_DIAG diag, ATL_CINT M, ATL_CINT N, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, TYPE *B, ATL_CINT ldb)
{
   ATL_ttrsm_amm_t pd;
   int P = (ATL_TP_PTR) ? ATL_TP_PTR->nthr : ATL_NTHREADS;
   int i;
   size_t sz, szA, szL;
   void *vp;

   P = Mjoin(PATL,tGetTrsmInfo)(&pd, P, TA, M, N, alpha);
   if (P < 2)
   {
      Mjoin(PATL,trsm)(AtlasLeft, uplo, TA, diag, M, N, alpha, A, lda, B, ldb);
      return(0);
   }
// printf("M=%d, N=%d, P=%d, mb=%d, nb=%d\n", M, N, P, pd.mb, pd.nb);
   pd.nablks = sz = ((pd.nmblks-1)*pd.nmblks)>>1;
   if (sz != pd.nablks)    /* if # of blocks overflow ints */
      return(1);              /* tell caller to recur until it doesn't */
   if (((size_t)pd.nmblks)*pd.nnblks != pd.nxblks)
      return(1);
   pd.TA = TA;
   pd.uplo = uplo;
   pd.diag = diag;
   pd.M = M;
   pd.N = N;
   pd.alpha = alpha;
   pd.lda = lda;
   pd.ldx = ldb;
   pd.A = A;
   pd.X = B;
   pd.blkszB = pd.mb * pd.nb;
   pd.blkszA = pd.mb * pd.mb;
   pd.panszC = (pd.nmblks-1) * pd.blkszB;
   szA = pd.nablks * pd.blkszA;
   pd.wsL = szL = pd.blkszB + pd.panszC;
   sz = szA + P*szL + pd.mu*pd.nu*pd.ku;
   sz = ATL_MulBySize(sz) + ATL_Cachelen;
   if ((sz>>20) > ATL_PTMAXMALLOC_MB)
      return(2);
   vp = malloc(sz);
   if (!vp)
      return(2);
   pd.AcpyDone = (pd.mb >= M);
   pd.wA = ATL_AlignPtr(vp);
   pd.w = pd.wA + szA;
/*
 * Select the number of cores to perform A copy: to many and they just
 * fight for the bus, to few and algorithm starts slower
 */
   if (P >= 8)
   {
      i = P>>2;
      i = (i > 4) ? 4 : i;
   }
   else if (P >= 4)
      i = P>>1;
   else
      i = P;
   pd.AblkCtr = ATL_SetAtomicCount(pd.nablks);  /* w/o glob, is in-order */
   pd.rhsCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nnblks,P), pd.nnblks, 0);
   pd.AcpyBV = ATL_NewBV(pd.nablks);
   pd.Acpymut = ATL_mutex_init();
/*   #define DEBUG1  */
   #ifdef DEBUG1
   {
      ATL_tpool_t *pp=ATL_TP_PTR;
      if (!pp)
         pp = ATL_NewThreadPool(1, 0, NULL);
      ATL_assert(pp);
      pp->PD = &pd;
      Mjoin(PATL,DoTrsm_amm)(pp, 0, 0);
      if (pp != ATL_TP_PTR)
         ATL_FreeThreadPool(pp);
   }
   #else
      ATL_goParallel(P, Mjoin(PATL,DoTrsm_amm), NULL, &pd, NULL);
   #endif
   ATL_FreeAtomicCount(pd.AblkCtr);
   ATL_FreeGlobalAtomicCount(pd.rhsCtr);
   ATL_FreeBV(pd.AcpyBV);
   ATL_mutex_free(pd.Acpymut);
   free(vp);
   return(0);
}

/*
 * Simple recursion for Left, Lower, Notrans case (used in LU)
 */
#include "atlas_tcacheedge.h"
#if CacheEdge > ATL_PTMAXMALLOC_MB*1024*1024
   #undef CacheEdge
#endif
#ifndef CacheEdge
   #define CacheEdge 524288
#endif

static void trsmREC_LLN
   (const int STOP_EARLY, const enum ATLAS_DIAG diag, ATL_CSZT M, ATL_CSZT N,
    const SCALAR alpha, const TYPE *A, ATL_CSZT lda, TYPE *B, ATL_CSZT ldb)
{
   ATL_CSZT Ml=M>>1, Mr = M-Ml;
   TYPE *B1=B+(Ml SHIFT);
   const TYPE *A10=A+(Ml SHIFT), *L11=A10+Ml*(lda SHIFT);
   #ifdef TCPLX
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
      TYPE NONE[2] = {ATL_rnone, ATL_rzero};
   #else
      #define ONE ATL_rone
      #define NONE ATL_rnone
   #endif
   size_t sz;  /* size(A) + 2*pan(B) */
/*
 * Stop recursion when A & panel of B fit in cache
 */
   sz = ((M*M)>>1) + M*(ATL_AMM_98KB<<1);
   sz = ATL_MulBySize(sz);
   if (STOP_EARLY || sz < CacheEdge)
   {
      if (!ttrsm_ammL(AtlasLower, AtlasNoTrans, diag, M, N, alpha,
                      A, lda, B, ldb))
         return;
/*
 *    If we can't allocate any space, call serial & hope it can work
 */
      if (M < 80)
      {
         Mjoin(PATL,trsm)(AtlasLeft, AtlasLower, AtlasNoTrans, diag, M, N,
                          alpha, A, lda, B, ldb);
         return;
      }
   }
   trsmREC_LLN(STOP_EARLY, diag, Ml, N, alpha, A, lda, B, ldb);
   Mjoin(PATL,tgemm)(AtlasNoTrans, AtlasNoTrans, Mr, N, Ml, NONE, A10, lda,
                     B, ldb, alpha, B1, ldb);
   trsmREC_LLN(STOP_EARLY, diag, Mr, N, ONE, L11, lda, B1, ldb);
}

static void trsmREC_LUT
   (const int STOP_EARLY, const enum ATLAS_DIAG diag, ATL_CSZT M, ATL_CSZT N,
    const SCALAR alpha, const TYPE *A, ATL_CSZT lda, TYPE *B, ATL_CSZT ldb)
{
   ATL_CSZT Ml=M>>1, Mr = M-Ml;
   TYPE *B1=B+(Ml SHIFT);
   const TYPE *A10=A+(Ml*lda SHIFT), *L11=A10+(Ml SHIFT);
   #ifdef TCPLX
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
      TYPE NONE[2] = {ATL_rnone, ATL_rzero};
   #else
      #define ONE ATL_rone
      #define NONE ATL_rnone
   #endif
   size_t sz;  /* size(A) + 2*pan(B) */
/*
 * Stop recursion when A & panel of B fit in cache
 */
   sz = ((M*M)>>1) + M*(ATL_AMM_98KB<<1);
   sz = ATL_MulBySize(sz);
   if (STOP_EARLY || sz < CacheEdge)
   {
      if (!ttrsm_ammL(AtlasUpper, AtlasTrans, diag, M, N, alpha,
                      A, lda, B, ldb))
         return;
/*
 *    If we can't allocate any space, call serial & hope it can work
 */
      if (M < 80)
      {
         Mjoin(PATL,trsm)(AtlasLeft, AtlasUpper, AtlasTrans, diag, M, N,
                          alpha, A, lda, B, ldb);
         return;
      }
   }
   trsmREC_LUT(STOP_EARLY, diag, Ml, N, alpha, A, lda, B, ldb);
   Mjoin(PATL,tgemm)(AtlasTrans, AtlasNoTrans, Mr, N, Ml, NONE, A10, lda,
                     B, ldb, alpha, B1, ldb);
   trsmREC_LUT(STOP_EARLY, diag, Mr, N, ONE, L11, lda, B1, ldb);
}

#ifndef TCPLX
   #undef ONE
   #undef NONE
#endif

#ifdef ATL_ARCH_XeonPHI
  #define ATL_SE 1
#else
  #define ATL_SE 0
#endif
int Mjoin(PATL,ttrsm_amm)
   (const enum ATLAS_SIDE side, const enum ATLAS_UPLO uplo,
    const enum ATLAS_TRANS TA, const enum ATLAS_DIAG diag,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, TYPE *B, ATL_CINT ldb)
{
   if (side == AtlasLeft)
   {
      if (uplo == AtlasLower && TA == AtlasNoTrans)
      {
         trsmREC_LLN(ATL_SE, diag, M, N, alpha, A, lda, B, ldb);
         return(0);
      }
      if (uplo == AtlasUpper && TA == AtlasTrans)
      {
         trsmREC_LUT(ATL_SE, diag, M, N, alpha, A, lda, B, ldb);
         return(0);
      }
   }
   return(1);
}
