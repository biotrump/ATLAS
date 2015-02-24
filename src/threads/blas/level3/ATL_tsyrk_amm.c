#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"
/*
 * Translates C coordinate in lower triangular matrix to Cblk #.
 * nm is the number of Mblocks,
 * (i,j) is the coordinate of the block
 */
#define Mcoord2cblk(i_, j_, nm_) ((((nm_)+(nm_)-j-1)*j)>>1 - (j_) + (i_) - 1)
/*
 * Translates block number b_ to (i,j) coordinates assuming first panel
 * has nm_ blocks in it (including diagonal block)
 */
#define Mcblk2coord(NM_, B_, I_, J_) \
{ \
   register int n_ = (NM_)-1, b_=(B_), j_; \
   for (j_=0; b_ >= n_; j_++) \
   { \
      b_ -= n_; \
      n_--; \
   } \
   (J_) = j_; \
   (I_) = j_ + b_ + 1; \
}

void Mjoin(PATL,CombSyrk_ammK)
(
   void *vp,          /* void ptr to ATL_GEMV_t struct given to threads */
   const int myrank,  /* my processor rank */
   const int hisrank  /* rank entry to be combined into mine */
)
{
   ATL_tsyrk_ammK_t *pd=vp;
   TYPE *myC, *hisC;
   const unsigned int incW = pd->kb*(pd->mb+pd->nb), N = pd->nb * pd->mb;

   myC = pd->w + pd->wsz*myrank + incW;
   hisC = pd->w + pd->wsz*hisrank + incW;
   Mjoin(PATL,axpy)(N, ATL_rone, hisC, 1, myC, 1);
}

static void DoSyrkK(unsigned int rank, ATL_tsyrk_ammK_t *pd)
{
   TYPE *wA, *wB, *wC;
   ammkern_t amm = pd->amm_b0;
   const cm2am_t a2blk=pd->a2blk, b2blk=pd->b2blk;
   const unsigned int mb=pd->mb, nb=pd->nb, kb=pd->kb, nkblks=pd->nkblks,
      N=pd->N, nmu=pd->nmu, nnu=pd->nnu;
   const size_t mulA = (pd->TA) ? 1 : pd->lda;
   int lda=pd->lda;
   int kctr;
   wA = pd->w + pd->wsz*rank;
   wB = wA + mb*kb;
   wC = wB + nb*kb;

   while ((kctr = ATL_DecGlobalAtomicCount(pd->KbCtr, rank)))
   {
      const int kblk = nkblks - kctr;
      const TYPE *a = pd->A + ((size_t)kblk)*mulA*kb;
      if (kblk != nkblks-1)  /* normal full-kb operation */
      {
         a2blk(kb, N, ATL_rone, a, lda, wA);
         b2blk(kb, N, ATL_rone, a, lda, wB);
         amm(nmu, nnu, kb, wA, wB, wC, wA, wB, wC);
      }
      else /* last block of size kb0 */
      {
         a2blk(pd->kb0, N, ATL_rone, a, lda, wA);
         b2blk(pd->kb0, N, ATL_rone, a, lda, wB);
         if (amm == pd->amm_b0)
            pd->ammK_b0(nmu, nnu, pd->KB0, wA, wB, wC, wA, wB, wC);
         else
            pd->ammK_b1(nmu, nnu, pd->KB0, wA, wB, wC, wA, wB, wC);
      }
      amm = pd->amm_b1;
   }
/*
 * If I did no work, zero my wrkspace so I don't screw up combine!
 */
   if (amm != pd->amm_b1)
      Mjoin(PATL,zero)(mb*nb, wC, 1);
}

void Mjoin(PATL,DoWork_syrk_amm_K)(ATL_LAUNCHSTRUCT_t *lp, void *vp)
{
   ATL_thread_t *tp = vp;
   DoSyrkK(tp->rank, lp->opstruct);
}
int Mjoin(PATL,tsyrk_amm_K)
(
   const enum ATLAS_UPLO Uplo,
   const enum ATLAS_TRANS Trans,
   ATL_CINT N,
   ATL_CINT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CINT lda,
   const SCALAR beta,
   TYPE *C,
   ATL_CINT ldc
)
{
   amminfo_t mminfo;
   ATL_tsyrk_ammK_t pd;
   ablk2cmat_t Mjoin(PATL,tGetSyammInfo_K)
      (amminfo_t *out, const int P, enum ATLAS_TRANS TA, ATL_CSZT N,ATL_CSZT K);
   int kb=ATL_AMM_MAXKB, nkb = K / ATL_AMM_MAXKB, P = ATL_NTHREADS;
   int ku, kr, mb, nb, mu, nu;
   size_t sz;
   void *vp=NULL;

   if (nkb < P)
   {
      kb = ATL_AMM_98KB;
      nkb = K / ATL_AMM_98KB;
      if (nkb < P)
      {
         nkb = K / ATL_AMM_66KB;
         kb = ATL_AMM_66KB;
      }
   }
   if (nkb < P)
   {
      if (nkb < 2)
      {
          Mjoin(PATL,syrk)(Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
          return(0);
      }
      P = nkb;
   }
   pd.blk2c_b0 = Mjoin(PATL,tGetSyammInfo_K)(&mminfo, P, Trans, N, kb);
   kb = mminfo.kb;
   nkb = K / kb;

   mu = mminfo.mu;
   nu = mminfo.nu;
   pd.nmu = (N+mu-1) / mu;
   pd.nnu = (N+nu-1) / nu;
   pd.mb = mb = pd.nmu*mu;
   pd.nb = nb = pd.nnu*nu;
   pd.kb = mminfo.kb;
   sz = ((((size_t)mb)*nb)<<1) + (mb+nb)*kb;
   pd.wsz = sz;
   sz = ATL_MulBySize(sz)*P;
   vp = malloc(sz+ATL_Cachelen);
   if (!vp)
      return(1);
   pd.w = ATL_AlignPtr(vp);
   kr = K - nkb*kb;
   pd.kb0 = pd.KB0 = kr;
   ku = mminfo.ku;
   if (!kr)
   {
      pd.kb0 = pd.KB0 = kb;
      pd.ammK_b0 = mminfo.amm_b0;
      pd.ammK_b1 = mminfo.amm_b1;
   }
   else
   {
      #if ATL_AMM_MAXKMAJ > 1
         if (ATL_AMMFLG_KMAJOR(mminfo.flag))
         {
            pd.KB0 = ((kr+ku-1)/ku)*ku;
            if (ATL_AMMFLG_KRUNTIME(mminfo.flag))
            {
               pd.ammK_b0 = mminfo.amm_b0;
               pd.ammK_b1 = mminfo.amm_b1;
            }
            else
            {
               pd.ammK_b0 = mminfo.amm_k1_b0;
               pd.ammK_b1 = mminfo.amm_k1_b1;
            }
         }
         else
      #endif
      {
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag) && kr == (kr/ku)*ku &&
             kr > mminfo.kbmin)
         {
               pd.ammK_b0 = mminfo.amm_b0;
               pd.ammK_b1 = mminfo.amm_b1;
         }
         else
         {
               pd.ammK_b0 = mminfo.amm_k1_b0;
               pd.ammK_b1 = mminfo.amm_k1_b1;
         }
      }
   }
   pd.amm_b0 = mminfo.amm_b0;
   pd.amm_b1 = mminfo.amm_b1;
   pd.blk2c_b1 = mminfo.Cblk2cm;
   pd.a2blk = mminfo.a2blk;
   pd.b2blk = mminfo.b2blk;
   pd.A = A;
   pd.C = C;
   pd.alpha = &alpha;
   pd.beta = &beta;
   pd.nkblks = (kr) ? nkb+1 : nkb;
   pd.KbCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nkblks, P), pd.nkblks, 0);
   pd.Cmut = ATL_mutex_init();
   pd.BETA_APPLIED = SCALAR_IS_ONE(beta);
   pd.LOWER = (Uplo == AtlasLower);
   pd.TA = (Trans == AtlasTrans);
   pd.N = N;
   pd.lda = lda;
   pd.ldc = ldc;

//   #define DEBUG1 1
   #ifdef DEBUG1
   {
      ATL_LAUNCHSTRUCT_t ls;
      ATL_thread_t ts;
      ts.rank = 0;
      ts.P = 1;
      ls.opstruct = &pd;
      Mjoin(PATL,DoWork_syrk_amm_K)(&ls, &ts);
   }
   #else
      ATL_goparallel(P, Mjoin(PATL,DoWork_syrk_amm_K), &pd,
                     Mjoin(PATL,CombSyrk_ammK));
   #endif
/*
 * Answer is written back to rank0's workspace, extract it & write to C
 */
   {
      TYPE *wC = pd.w+kb*(mb+nb), *w = wC + mb*nb, *c = C;
/*
 *    Put it into block-major storage in w
 */
      pd.blk2c_b0(N, N, ATL_rone, wC, ATL_rzero, w, N);
/*
 *    Now copy out only upper or lower portion
 */
      if (pd.LOWER)
      {
         int j;
         for (j=0; j < N; j++, c += ldc+1, w += N+1)
             Mjoin(PATL,axpby)(N-j, alpha, w, 1, beta, c, 1);
      }
      else
      {
         int j;
         for (j=0; j < N; j++, c += ldc, w += N)
             Mjoin(PATL,axpby)(j+1, alpha, w, 1, beta, c, 1);
      }
   }
   free(vp);
   ATL_mutex_free(pd.Cmut);
   ATL_FreeGlobalAtomicCount(pd.KbCtr);
   return(0);
}

/*
 * On PHI, we deal out 4 Cblks at once, so that contexts on same core share
 * either A (Upper) or B (Lower)
 */
#ifdef ATL_PHI_SORTED
   #define DEBUG2 1
/*
 * RETURNS: >=0 : block to work on, -1: no work left, -2, no work avail
 */
static int FindContextChunkC(ATL_tsyrk_ammN_t *pd)
{
/*
 * While we haven't finished copying all of A, must check dependencies
 * before giving out work
 */
   if (ATL_FindFirstUnsetBitBV(pd->cpydonBV, 0) != -1)
   {
      const unsigned int ndiag=pd->ndiag, ncblks=pd->ncblks, LOWER=pd->LOWER,
                         ncntxts=pd->ncntxts;
      int ib = -1, bb=0, i, j, k, mul;

      do
      {
         ib = ATL_FindFirstUnsetBitBV(pd->cblkBV, bb);
         if (ib == -1)
            return((bb) ? -2 : -1);
         Mcblk2coord(ndiag, ib, i, j);
         k = i-j-1;

         if (ncntxts != 3)
            mul = !(k&(ncntxts-1));
         else
            mul = (k == (k/3)*3);
         if (mul && ATL_IsBitSetBV(pd->cpydonBV, i) &&
             ATL_IsBitSetBV(pd->cpydonBV, j))
         {
            int  nl = ndiag-i;
            nl = Mmin(nl, ncntxts);
/*
 *          If only one row/col left in dir of context, then OK
 */
            if (nl == 1)
               return(ib);
/*
 *          If >1 left, must check of 2nd context's operand ready as well
 */
            else
            {
               if (ATL_IsBitSetBV(pd->cpydonBV, i+1))
               {
                  if (nl == 2)
                     return(ib);
                  if (ATL_IsBitSetBV(pd->cpydonBV, i+2))
                  {
                     if (nl == 3)
                        return(ib);
                     if (ATL_IsBitSetBV(pd->cpydonBV, i+3))
                        return(ib);
                  }
               }
            }
         }     /* end ifs checking if A/B blocks have been copied */
         bb = ib + 1;
      }
      while (bb < ncblks);
      return(-2);
   }  /* end if I still am copying A */
/*
 * If everything has been copied, then just take the first available
 * chunk.  Since we have dealt out all other blocks in chunks, we know
 * any remaining chunks are full or are partial due to N
 */
   return(ATL_FindFirstUnsetBitBV(pd->cblkBV, 0));
}

static int GetCBlkChunk(const int rank, ATL_tsyrk_ammN_t *pd)
{
   int ib, bb=0, i, j, k;
   const unsigned int ndiag=pd->ndiag, ncblks=pd->ncblks, ncntxts=pd->ncntxts;
/*
 *  We will exit loop when we run out of panels (ib=-1), or we get a valid
 *  block (ib >= 0).  If we have a valid block, we'll also have cwmut locked.
 */
   while((ib = FindContextChunkC(pd)) != -1)
   {
      if (ib == -2)
      {
         ATL_thread_yield();
         continue;
      }
/*
 *    OK, there's a potential panel, so lock the data structure and if there
 *    is still a panel to be had, take it
 */
      ATL_mutex_lock(pd->cwmut);
      ib = FindContextChunkC(pd);
      if (ib >= 0)
         break;
      ATL_mutex_unlock(pd->cwmut);
   }   /* exit loop wt valid ib & cwmut
/*
 * If we are out of C blocks, inform thread (cwmut is not locked!)
 */
   if (ib == -1)
      return(-1);
/*
 * Otherwise, we exited loop with ib>=0 and the cwmut locked
 * Lower deals out blks from column panel, Upper from row panel;
 * Get 4 of them in a row, unless we reach end of panel; in this case
 * 1-3 contexts will be inactive during the computation of the cblk
 * to avoid pollution of the cache
 */
   Mcblk2coord(ndiag, ib, i, j);
   #ifdef DEBUG2
      ATL_assert(!ATL_IsBitSetBV(pd->cblkBV,ib));
   #endif
   ATL_SetBitBV(pd->cblkBV, ib);
   if (i+1 < ndiag)
   {
      #ifdef DEBUG2
         ATL_assert(!ATL_IsBitSetBV(pd->cblkBV,ib+1));
      #endif
      ATL_SetBitBV(pd->cblkBV, ib+1);
      if (i+2 < ndiag && ncntxts > 2)
      {
         #ifdef DEBUG2
            ATL_assert(!ATL_IsBitSetBV(pd->cblkBV,ib+2));
         #endif
         ATL_SetBitBV(pd->cblkBV, ib+2);
         if (i+3 < ndiag && ncntxts > 3)
         {
            #ifdef DEBUG2
               ATL_assert(!ATL_IsBitSetBV(pd->cblkBV,ib+3));
            #endif
            ATL_SetBitBV(pd->cblkBV, ib+3);
         }
      }
   }
   ATL_mutex_unlock(pd->cwmut);
   return(ib);
}
#endif

static int GetCBlk(const int rank, ATL_tsyrk_ammN_t *pd)
{
   int ib, bb=0;
   const unsigned int ndiag=pd->ndiag, ncblks=pd->ncblks;
   do
   {
      int i, j, k;
/*
 *    If all copying has finished, we just take the first available block
 */
      if (ATL_FindFirstUnsetBitBV(pd->cpydonBV, 0) == -1)
         pd->cpydone = 1;
      if (pd->cpydone)
      {
         bb = ATL_FindFirstUnsetBitBV(pd->cblkBV, 0);
         if (bb == -1)
            return(-1);
         ATL_mutex_lock(pd->cwmut);
         ib = ATL_FindFirstUnsetBitBV(pd->cblkBV, bb);
         if (ib != -1)
            ATL_SetBitBV(pd->cblkBV, ib);
         ATL_mutex_unlock(pd->cwmut);   /* unlock mutex protecting cblkBV */
         return(ib);
      }
      ib = ATL_FindFirstUnsetBitBV(pd->cblkBV, bb);
      if (ib == -1)
      {
         if (!bb)
            return(-1);
         else
            bb = 0;
         continue;
      }
      bb = ib+1;          /* don't check this again until we cycle */
      if (bb == ncblks)   /* if at end of bitvec */
         bb = 0;          /* start looking again at beginning */
/*
 *    There's an unassigned block, check if it's A & A' blks have been copied
 */
      Mcblk2coord(ndiag, ib, i, j);
//      fprintf(stderr, "%d:ib=%d, i=%d, j=%d\n", rank, ib, i, j);
      if (!ATL_IsBitSetBV(pd->cpydonBV, i))
         continue;
      if (!ATL_IsBitSetBV(pd->cpydonBV, j))
         continue;
/*
 *    Dep look good, now seize mutex and take it for real
 */
      ATL_mutex_lock(pd->cwmut);
      k = ATL_FindFirstUnsetBitBV(pd->cblkBV, ib);
/*
 *    If nobody else got the proven-good block before me, take it
 */
      if (k == ib)
      {
         ATL_SetBitBV(pd->cblkBV, ib);  /* I take it */
         ATL_mutex_unlock(pd->cwmut);   /* unlock mutex protecting cblkBV */
         return(ib);                    /* caller can compute this block */
      }
      ATL_mutex_unlock(pd->cwmut);
   }
   while(1);
   return(-1);
}

/*
 * computes (i,j) non-diagonal block of C
 */
#ifndef ATL_PHI_SORTED2
static void DoCblk(const int rank, ATL_tsyrk_ammN_t *pd, TYPE *wC, int i, int j)
#else
static void DoCblk(const int rank, ATL_tsyrk_ammN_t *pd, TYPE *wC, int i, int j,
                   cont int crank)
#endif
{
   const ammkern_t amm = pd->amm_b1;
   const unsigned int nkblks=pd->nkblks, bs=pd->blkszA, kb=pd->kb, NB=pd->nb;
   unsigned int nmu, nnu, mb, nb;
   #ifdef ATL_PHI_SORTED2
      unsigned int nnuC, nr;  /* nnu / ncntxts */
   #endif
   const TYPE *wA, *wB, *wAn, *wBn;
   TYPE *c;
   int k;

   if (!(pd->LOWER))
   {
      k = i;
      i = j;
      j = k;
   }
   if (j != pd->ndiag-1)
   {
      nnu = pd->nnu;
      nb = pd->nb;
   }
   else
   {
      nnu = pd->nnuf;
      nb = pd->nbf;
   }
   if (i != pd->ndiag-1)
   {
      nmu = pd->nmu;
      mb = pd->nb;
   }
   else
   {
      nmu = pd->nmuf;
      mb = pd->nbf;
   }
   wA = pd->wA + i*pd->panszA;
   wB = pd->wAt + j*pd->panszA;
   wA = pd->wA + i*pd->panszA;
   wB = pd->wAt + j*pd->panszA;
   #ifdef ATL_PHI_SORTED2
      nnuC = nnu / pd->ncntxts;
      nr = nnu - nnuC * pd->ncntxts;
      if (crank == 1)
         pB +=
      else if (crank == 2)
      else if (crank == 3)
      pB += bs / pd->ncntxts;
      if (crank < nr)
         nnuC++;
   #endif
   wAn = wA + bs;
   wBn = wB + bs;
   #ifdef DEBUG2
      if (!ATL_IsBitSetBV(pd->cpydonBV, i) || !ATL_IsBitSetBV(pd->cpydonBV, j))
          fprintf(stderr, "%d: ndiag=%d, i=%d, j=%d\n", rank, pd->ndiag, i, j);
      ATL_assert(ATL_IsBitSetBV(pd->cpydonBV, i));
      ATL_assert(ATL_IsBitSetBV(pd->cpydonBV, j));
   #endif
   pd->ammK(nmu, nnu, pd->KB0, wA, wB, wC, wAn, wBn, wC);
   for (k=1; k < nkblks; k++)
   {
      wA = wAn;
      wB = wBn;
      wAn += bs;
      wBn += bs;
      amm(nmu, nnu, kb, wA, wB, wC, wAn, wBn, wC);
   }
   pd->blk2c(mb, nb, *(pd->alpha), wC, *(pd->beta),
             pd->C+ NB*(j*(size_t)(pd->ldc) + i), pd->ldc);
}

#ifdef ATL_PHI_SORTED
static void DoNonDiagChunk(const int rank, ATL_tsyrk_ammN_t *pd)
{
   const unsigned int ndiag = pd->ndiag, LOWER=pd->LOWER, ncntxts=pd->ncntxts;
   const unsigned crank = rank / pd->ncores;  /* my context rank */
   int ic;

   TYPE *wC = pd->wC + (rank+rank)*(pd->nbnb);
   volatile int *chkin;

   ic = rank - crank*pd->ncores;  /* What core am I on? */
   chkin = (volatile int*) (((volatile char *)pd->chkin)+ATL_Cachelen*ic);
/*
 * Each core waits for all its contexts to show up before doing work, so that
 * we don't grab work well in advance of all contexts arriving
 * 0 makes sure everyone reads the initial chkin[0] before proceeding, since
 * the following sync switches who writes first
 */
   if (!crank)
   {
      if (ncntxts == 4)
      {
         while (chkin[1] != -2 || chkin[2] != -2 || chkin[3] != -2)
            ATL_thread_yield();
         chkin[4] = chkin[5] = chkin[6] = -1;
         *chkin = -2;
         while (chkin[4] != -8 || chkin[5] != -8 || chkin[6] != -8)
            ATL_thread_yield();
      }
      else if (ncntxts == 3)
      {
         while (chkin[1] != -2 || chkin[2] != -2)
            ATL_thread_yield();
         chkin[4] = chkin[5] = -1;
         *chkin = -2;
         while (chkin[4] != -8 || chkin[5] != -8)
            ATL_thread_yield();
      }
      else /* if (ncntxts == 2) */
      {
         while (chkin[1] != -2)
            ATL_thread_yield();
         chkin[4] = -1;
         *chkin = -2;
         while (chkin[4] != -8)
            ATL_thread_yield();
      }
   }
   else
   {
      ic = chkin[crank] = -2;
      while (chkin[0] != -2)
         ATL_thread_yield();
      chkin[crank+3] = -8;
   }

//   fprintf(stderr, "%d: START DoNonDiag: ndiag=%d, ncblks=%d\n",
//           rank, ndiag, pd->ncblks);
   while (1)
   {
      int i, j;

      if (!crank)
      {
         ic = GetCBlkChunk(rank, pd);
//         fprintf(stderr, "\n%d: NEW ic=%d\n", rank, ic);
         *chkin = ic;
         if (ncntxts == 4)
         {
            while (chkin[1] != ic || chkin[2] != ic || chkin[3] != ic)
               ATL_thread_yield();
         }
         else if (ncntxts == 3)
         {
            while (chkin[1] != ic || chkin[2] != ic)
               ATL_thread_yield();
         }
         else /* if (ncntxts == 2) */
         {
            while (chkin[1] != ic)
               ATL_thread_yield();
         }
      }
      else
      {
         while (chkin[0] == ic)
            ATL_thread_yield();
         chkin[crank] = ic = chkin[0];
      }
      if (ic == -1)
         return;
      Mcblk2coord(ndiag, ic, i, j);
      #ifdef DEBUG2
         ATL_assert(ATL_IsBitSetBV(pd->cblkBV, ic));
      #endif
      i += crank;
      if (i < ndiag)
         DoCblk(rank, pd, wC, i, j);
   }
}
#endif

/*
 * For non-diagonal work, we count the number of non-diagonal blocks of C,
 * which is initially stored in the ncblks variable, which is protected
 * the cwmut mutex, which also protects the cblkBV, which is a ncblk-len BV.
 * A unset bit means that particular non-diagonal block has not yet been
 * assigned to a thread, while a 1 means it has.  The mutex also protects
 * the cpydone variable, which is set to 1 when cpydonBV has all bits set.
 * So, threads wanting to do non-diagonal work will find the first unset
 * bit in cblkBV, and then translate that to a (i,j) C block coordinate.
 * If the ith & jth bits are both set in cpydonBV (or cpydone is set), then
 * they will take that block as their own to do (in this phase, each thread
 * gets an individual block of C to do) by setting the bit in cblkBV.
 * When all bits are set in cblkBV, then all work has been dealt out, and
 * threads will exit once they say there is no more work to do.
 * The master process joins all created threads, and can delete data structures
 * safely after all joins succeed.
 */
#ifndef ATL_PHI_SORTED2
static void DoNonDiag(const int rank, ATL_tsyrk_ammN_t *pd)
{
   const unsigned int ndiag = pd->ndiag;
   int ic;
   TYPE *wC = pd->wC + (rank+rank)*(pd->nbnb);

//   fprintf(stderr, "%d: START DoNonDiag: ndiag=%d, ncblks=%d\n",
//           rank, ndiag, pd->ncblks);
   while ((ic = GetCBlk(rank, pd)) != -1)
   {
      int i, j;
      Mcblk2coord(ndiag, ic, i, j);
//      fprintf(stderr, "%d: ic=%d, i=%d, j=%d\n", rank, ic, i, j);
      DoCblk(rank, pd, wC, i, j);
   }
}
#else
static void DoNonDiag(const int rank, ATL_tsyrk_ammN_t *pd)
{
   const unsigned int ndiag = pd->ndiag, crank = rank/pd->ncores;
   const unsigned int ncntxts=pd->ncntxts;
   unsigned int nnuC, nr;
   int ic, k;
   TYPE *wC = pd->wC + (rank+rank)*(pd->nbnb);
   volatile int *chkin;
   int mycnt;
   chkin = (volatile int*) (((volatile char *)pd->chkin)+ATL_Cachelen*ic);
   mycnt = chkin[crank];
   nnuC = pd->nnu / pd->ncntxts;
   nr = pd->nnu - nnuC*pd->ncntxts;
   for (k=ic=0; ic < crank; ic++)
      k += (ic < nr) ? nnuC+1 : nnuC;
   wC += k*mu*nb;   /* k*mu*nu*(nb/nu) = k*mu*nb */
   // HERE HERE HERE

//   fprintf(stderr, "%d: START DoNonDiag: ndiag=%d, ncblks=%d\n",
//           rank, ndiag, pd->ncblks);
   while (1)
   {
      int i, j;
      if (crank)  /* I'm a worker context */
      {
         while (*chkin == mycnt)
            ATL_thread_yield();
         ic = mycnt = chkin[crank] = *chkin;
      }
      else /* I'm master context querying for block to work on */
      {
         ic = GetCBlk(rank, pd);
         if (ncntxts == 4)
         {
            while (chkin[1] != mycnt || chkin[2] != mycnt || chkin[3] != mycnt)
               ATL_thread_yield();
            mycnt = chkin[0] = ic;
         }
      }
      if (ic == -1)   /* if out of work */
         break;       /* can quit function */
      Mcblk2coord(ndiag, ic, i, j);
      DoCblk(rank, pd, wC, i, j);
   }
}
#endif

static void DoBlkWtCpy
(
   unsigned const int rank,   /* my thread rank */
   ATL_tsyrk_ammN_t *pd,      /* problem def structure */
   unsigned const int dblk,   /* diagonal block of C to compute */
   unsigned const int kblk,   /* kblk to compute */
   unsigned const int b,      /* actual size of block w/o padding */
   unsigned const int nmu,    /* m/mu for this C blk */
   unsigned const int nnu,    /* n/nu for this C blk */
   ammkern_t amm,             /* kern to use for all non-kb0 blocks */
   TYPE *wC                   /* my private workspace */
)
{
   const int nb=pd->nb, kb=pd->kb;
   const TYPE *A;             /* ptr to original matrix */
   TYPE *wA, *wB;             /* ptr to copy space */
   size_t i, k;
   int kk;
// fprintf(stderr, "%d: dblk=%d, kblk=%d\n", rank, dblk, kblk);
   i = dblk * nb;
   k = (kblk) ? pd->kb0 + (kblk-1) * kb : 0;
   if (pd->TA)  /* A is transposed */
   {
      A = pd->A + k + i*pd->lda;
   }
   else        /* A is Notranspose */
   {
      A = pd->A + k*pd->lda + i;
   }
   wA = pd->wA + dblk*(pd->panszA) + kblk*(pd->blkszA);
   wB = pd->wAt + dblk*(pd->panszA) + kblk*(pd->blkszA);
   if (kblk)
   {
      pd->a2blk(kb, b, ATL_rone, A, pd->lda, wA);
      pd->b2blk(kb, b, ATL_rone, A, pd->lda, wB);
/*
 *    If there are any off-diag blks awaiting the copy, signal copy of this
 *    K-panel of A/A' is complete iff I just copied the last block
 */
      if (pd->ncblks)
      {
         if (ATL_DecGlobalAtomicCountDown((pd->KdonCtr)[dblk], rank) == 1)
         {
            ATL_mutex_lock(pd->cdmut);
            ATL_SetBitBV(pd->cpydonBV, dblk);
            ATL_mutex_unlock(pd->cdmut);
         }
      }
      amm(nmu, nnu, kb, wA, wB, wC, wA, wB, wC);
   }
   else
   {
      pd->a2blk(pd->kb0, b, ATL_rone, A, pd->lda, wA);
      pd->b2blk(pd->kb0, b, ATL_rone, A, pd->lda, wB);
/*
 *    If there are any off-diag blks awaiting the copy, signal copy of this
 *    K-panel of A/A' is complete iff I just copied the last block
 */
      if (pd->ncblks)
      {
         if (ATL_DecGlobalAtomicCountDown((pd->KdonCtr)[dblk], rank) == 1)
         {
            ATL_mutex_lock(pd->cdmut);
            ATL_SetBitBV(pd->cpydonBV, dblk);
            ATL_mutex_unlock(pd->cdmut);
         }
      }
      if (amm == pd->amm_b0)   /* ammK works for beta=0, so OK */
         pd->ammK(nmu, nnu, pd->KB0, wA, wB, wC, wA, wB, wC);
/*
 *    If using ammK in middle of operation, it handles only beta=0, so do
 *    operation in extra workspace and then add result into running sum
 */
      else
      {
         pd->ammK(nmu, nnu, pd->KB0, wA, wB, wC+pd->nbnb, wA, wB, wC+pd->nbnb);
         Mjoin(PATL,axpy)(pd->nbnb, ATL_rone, wC+pd->nbnb, 1, wC, 1);
      }
   }
}

static void DoBlksWtCpy
(
   unsigned const int rank,   /* my thread rank */
   ATL_tsyrk_ammN_t *pd,      /* problem def structure */
   unsigned const int dblk,   /* diagonal block of C to compute */
   unsigned int kctr,         /* non-zero Kbeg ctr */
   TYPE *wC                   /* my private workspace */
)
{
   const int TRANS = (pd->TA == AtlasTrans);
   int kblk = pd->nkblks - kctr;
   const int b = (dblk == pd->ndiag-1) ? pd->nbf : pd->nb;
   const int nmu = (dblk == pd->ndiag-1) ? pd->nmuf : pd->nmu;
   const int nnu = (dblk == pd->ndiag-1) ? pd->nnuf : pd->nnu;
   TYPE *wA, *wB, *w, *c;
   TYPE beta = ATL_rone;
   DoBlkWtCpy(rank, pd, dblk, kblk, b, nmu, nnu, pd->amm_b0, wC);
   while (kctr = ATL_DecGlobalAtomicCount((pd->KbegCtr)[dblk], rank))
   {
      kblk = pd->nkblks - kctr;
      DoBlkWtCpy(rank, pd, dblk, kblk, b, nmu, nnu, pd->amm_b1, wC);
   }
/*
 * Since I'm done with this blk of C, copy it to block-major storage, and
 * scale by alpha (blk2c won't access only triangle, so must cpy first)
 */
   w = wC + pd->nbnb;
   pd->blk2d(b, b, *(pd->alpha), wC, ATL_rzero, w, b);
/*
 * Now, seize mutex for diagonal block of original C, and copy back out
 * only above/below diagonal
 */
   c = pd->C + dblk*(pd->nb)*(pd->ldc+1);
   ATL_mutex_lock(pd->Cdmuts[dblk]);
   if (!ATL_IsBitSetBV(pd->dbetaBV, dblk)) /* if I apply beta */
   {
      ATL_SetBitBV(pd->dbetaBV, dblk); /* tell rest of thrs don't apply beta */
      beta = *(pd->beta);              /* because I'm going to do it */
   }
   if (pd->LOWER)
   {
      int j;
      for (j=0; j < b; j++, c += pd->ldc+1, w += b+1)
          Mjoin(PATL,axpby)(b-j, ATL_rone, w, 1, beta, c, 1);
   }
   else
   {
      int j;
      for (j=0; j < b; j++, c += pd->ldc, w += b)
          Mjoin(PATL,axpby)(j+1, ATL_rone, w, 1, beta, c, 1);
   }
   ATL_mutex_unlock(pd->Cdmuts[dblk]);
}

/*
 * In the first phase, we work only on diagonal blocks, while copying both
 * A & A'.  For diag work, we parallelize both N & K dims so that the copy
 * is done as quickly as possible.  Threads coming in first choose differing
 * diag blks; diagonal blocks are dealt out cheaply using the dCtr global
 * counter (which starts at nnblks == ndiag).
 * Once all diagonal blocks are dealt out, new threads will start using
 * the atomic ctr array KbegCtr array to share K work for each diagonal.
 * both KbegCtr & KdonCtr are nnblk-len arrays of atomic counters.  Each
 * counter starts at nkblks.  Once the block pointed to by KbegCtr is
 * completely copied, the copying array increments the KdonCtr.  Only one
 * core per diag will get KdonCtr == 0 after doing his copy, and this
 * core will seize cdmut mutex in order to set the appropriate bit in
 * cpydonBV, which is a nnblks-length bit vector.  If the kth bit is set,
 * that means the ith row of A & jth col of A' has been copied.
 * Once a thread gets KbegCtr for a particular diag of 0, it means there's
 * no more work for this block of C, and so it will seize the appropriate
 * Cdmuts mutex which protects each diagonal block of C, and write its
 * finished contribution out to C.  The first such thread to ever seize
 * the mutex will scope dbetaBV to find this diagonal block needs beta applied;
 * while later threads will use beta=1.
 * Eventually, all diagonal work is finished, and the first processor to
 * get 0 for all dCtr & KbegCtr requests will set NODWORK=1, so later
 * threads don't have to query all the counters to know they should proceed
 * to non-diagonal work.
 */
static void DoDiag(const int rank, ATL_tsyrk_ammN_t *pd)
{
   int DIAG=1, k;
   TYPE *myC = pd->wC + (rank+rank)*(pd->nbnb);
   while (!(pd->NODWORK))
   {
       int d=0;

/*
 *     Find which diagonal block to work on, and then which k blk to use
 */
       if (DIAG)
       {
          d = ATL_DecGlobalAtomicCount(pd->dCtr, rank);
          if (d)
          {
             k = ATL_DecGlobalAtomicCount((pd->KbegCtr)[pd->ndiag - d], rank);
             if (!k)     /* if no more K work to do */
                d = 0;   /* can't work on this diag after all */
          }
       }
/*
 *     If all diagonal blocks currently being worked on by threads, find
 *     one that I can help with.
 */
       if (!d)
       {
          unsigned int i, n=pd->ndiag;
          DIAG = 0;
          for (i=0; i < n; i++)
          {
             unsigned int j = (i+rank)%n;
             k = ATL_DecGlobalAtomicCount((pd->KbegCtr)[j], rank);
             d = n-j;
             if (k)
                goto FOUNDDK;
          }
          pd->NODWORK = 1;    /* no work left to assign */
          return;             /* so done with copy of A/A' & this func */
       }
/*
 *     If I reach here, I've got a valid d & k;  and I'll call a routine
 *     that continues to grab blocks from this diag along K until all K
 *     is done; it will then write the answer back to the original C, and
 *     return to this loop to see if it can help with another diag.
 */
       FOUNDDK:
          DoBlksWtCpy(rank, pd, pd->ndiag-d, k, myC);
   }
}

void Mjoin(PATL,DoWork_syrk_amm)(ATL_LAUNCHSTRUCT_t *lp, void *vp)
{
   ATL_thread_t *tp = vp;
   ATL_tsyrk_ammN_t *pd = lp->opstruct;
   if (!(pd->NODWORK))
      DoDiag(tp->rank, pd);
   if (pd->ncblks)
   {
      #ifdef ATL_PHI_SORTED
         if (pd->ncntxts > 1)
            DoNonDiagChunk(tp->rank, pd);
         else
      #endif
      if (ATL_FindFirstUnsetBitBV(pd->cblkBV, 0) != -1)
         DoNonDiag(tp->rank, pd);
   }
}
/*
 * SYRK where all parallelism comes from blocks of N, built atop amm directly
 *    if (TA == AtlasNoTrans)
 *       C = alpha * A*A' + beta*C
 *    else
 *       C = alpha * A'*A + beta*C
 *    C is an upper or lower symmetric NxN matrix,
 *    A is a dense rectangular NxK (NoTrans) or KxN (Trans) matrix
 * RETURNS: 0 if operation performed, non-zero otherwise.
 *   This routine assumes it can copy all of A up-front to simplify parallelism.
 *   Will return non-zero if memory cannot be allocated, on the assumption
 *   it is called from recursive implementation that can recur until malloc
 *   succeeds.
 */
int Mjoin(PATL,tsyrk_amm_N)
(
   const enum ATLAS_UPLO Uplo,
   const enum ATLAS_TRANS Trans,
   ATL_CINT N,
   ATL_CINT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CINT lda,
   const SCALAR beta,
   TYPE *C,
   ATL_CINT ldc
)
{
   int NB, nb, mu, nu, ku, nfdiag, ndiag, nkblks, ncnt, nkcnt, blkszA, nbf;
   int ncblks, kb, kb0, KB0, i, P;
   size_t sz, szA;
   ATL_tsyrk_ammN_t pd;
   amminfo_t mminfo;
   void *vp;
   ablk2cmat_t Mjoin(PATL,tGetSyammInfo)
      (amminfo_t *out, const int P, enum ATLAS_TRANS TA, ATL_CSZT N,
       ATL_CSZT K, const SCALAR alpha, const SCALAR beta);

   #ifdef ATL_AMM_98NB
      ncblks = N / ATL_AMM_98NB;
   #else
      ncblks = N / ATL_AMM_MAXNB;
   #endif
   ncblks = (ncblks*(ncblks-1))>>1;
   #ifdef ATL_PHI_SORTED2
      P = (ncblks >= ATL_NTHREADS>>2) ? ATL_NTHREADS : (ATL_NTHREADS>>1);
      pd.blk2d = Mjoin(PATL,tGetSyammInfo)(&mminfo, P, Trans, N, K, alpha,beta);
      nb = mminfo.nb;
      nu = mminfo.nu;
      pd.nnu = nb / nu;
      ncblks = N / nb;
      ncblks = (ncblks*(ncblks-1))>>1;
      if (ncblks >= ATL_NTHREADS>>2)
      {
         pd.ncores = ATL_NTHREADS>>2;
         pd.ncntxts = (pd.nnu >= 4) ? 4 : pd.nnu;
      }
      P = ncblks<<2;
      if (P
      if (ncblks <
      pd.ncores = ATL_NTHREADS>>2;
      if (P >= ATL_NTHREADS>>2)
         pd.ncntxts = 4;
      else if (P >= ATL_NTHREADS
   #elif defined(ATL_PHI_SORTED)
      if (ncblks >= ATL_NTHREADS)
      {
         P = ATL_NTHREADS;
         pd.ncores = ATL_NTHREADS>>2;
         pd.ncntxts = 4;
      }
      else if (ncblks+8 >= 3*(ATL_NTHREADS>>2))
      {
         pd.ncores = ATL_NTHREADS>>2;
         P = 3*pd.ncores;
         pd.ncntxts = 3;
      }
      else if (ncblks+4 >= (ATL_NTHREADS>>1))
      {
         P = ATL_NTHREADS>>1;
         pd.ncores = ATL_NTHREADS>>2;
         pd.ncntxts = 2;
      }
      else
      {
         pd.ncores = P = ncblks;
         pd.ncntxts = 1;
      }
   #else
      P = (ncblks < ATL_NTHREADS) ? ncblks : ATL_NTHREADS;
   #endif
   if (P < 2)
   {
      Mjoin(PATL,syrk)(Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
      return(0);
   }
   pd.blk2d = Mjoin(PATL,tGetSyammInfo)(&mminfo, P, Trans, N, K, alpha, beta);
   pd.blk2c = mminfo.Cblk2cm;
   pd.a2blk = mminfo.a2blk;
   pd.b2blk = mminfo.b2blk;
   nb = mminfo.nb;
   kb = mminfo.kb;
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   nfdiag = N/nb;
   nbf = N - nfdiag*nb;
   ndiag = (nbf) ? nfdiag+1 : nfdiag;
   ncblks = ((ndiag-1)*ndiag)>>1;
/*
 * Set up parallel data structure
 */
   nkblks = (K+kb-1)/kb;
   KB0 = kb0 = K - (K/kb)*kb;
   if (!kb0)
   {
      kb0 = KB0 = kb;
      pd.ammK = mminfo.amm_b0;
   }
   else
   {
      #if ATL_AMM_MAXKMAJ > 1
         if (ATL_AMMFLG_KMAJOR(mminfo.flag))
         {
            KB0 = ((kb0+ku-1)/ku)*ku;
            if (ATL_AMMFLG_KRUNTIME(mminfo.flag))
               pd.ammK = mminfo.amm_b0;
            else
               pd.ammK = mminfo.amm_k1_b0;
         }
         else
      #endif
      {
         pd.ammK = mminfo.amm_k1_b0;
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag) && kb0 == (kb0/ku)*ku &&
             kb0 > mminfo.kbmin)
            pd.ammK = mminfo.amm_b0;
      }
   }
   pd.amm_b0 = mminfo.amm_b0;
   pd.amm_b1 = mminfo.amm_b1;
   blkszA = kb * nb;
   szA = ndiag*nkblks*blkszA;                 /* space for A & A' */
   sz = (ndiag+ndiag+ndiag)*sizeof(void*);    /* space for Ctr & mut arrays */
   sz += ATL_MulBySize(nb)*nb*2*P; /* local C wrkspc */
   sz += ATL_MulBySize(szA+szA);              /* A/At workspace */
   sz += 3*ATL_Cachelen;                      /* room for alignment */
   #ifdef ATL_PHI_SORTED
      if (pd.ncntxts > 1)
         sz += ATL_Cachelen*(pd.ncores+1);
   #endif
   if ((sz>>20) > ATL_PTMAXMALLOC_MB)       /* if I'm over malloc limit */
      return(2);                            /* return and recur on K */
   vp = malloc(sz);
   if (!vp)
      return(2);
   #ifdef ATL_PHI_SORTED
   if (pd.ncntxts > 1)
   {
      pd.chkin = (volatile int*)ATL_AlignPtr(vp);
      pd.KbegCtr = (void*)(((char*)pd.chkin)+pd.ncores*ATL_Cachelen);
      for (i=0; i < pd.ncores; i++)
      {
         int *ip;
         char *cp = (char*)pd.chkin;

         cp += i*ATL_Cachelen;
         ip = (int*) cp;
         *ip = ip[1] = ip[2] = ip[3] = -128;
      }
   }
   else
   #endif
      pd.KbegCtr = vp;
   pd.alpha = &alpha;
   pd.beta = &beta;
   pd.TA = (Trans == AtlasTrans);
   pd.LOWER = (Uplo == AtlasLower);
   pd.ndiag = ndiag;
   pd.ncblks = ncblks;
   pd.blkszA = blkszA;
   pd.panszA = nkblks*blkszA;
   pd.nkblks = nkblks;
   pd.nb = nb;
   pd.nbnb = nb*nb;
   pd.kb = kb;
   pd.kb0 = kb0;
   pd.KB0 = KB0;
   pd.nnu = nb / nu;
   pd.nmu = nb / mu;
   if (nbf)
   {
      pd.nbf = nbf;
      pd.nmuf = (nbf+mu-1)/mu;
      pd.nnuf = (nbf+nu-1)/nu;
      pd.Mf = pd.nmuf*mu;
      pd.Nf = pd.nnuf*nu;
   }
   else
   {
      pd.nbf = pd.Mf = pd.Nf = nb;
      pd.nmuf = pd.nmu;
      pd.nnuf = pd.nnu;
   }

   pd.cpydone = pd.NODWORK = 0;
   pd.dCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(ndiag, P), ndiag, 0);
   nkcnt = ATL_EstNctr(nkblks, P);
   pd.KdonCtr = pd.KbegCtr + ndiag;
   pd.Cdmuts = pd.KdonCtr + ndiag;
   pd.wA = (TYPE*)(pd.Cdmuts+ndiag);
   pd.wA = ATL_AlignPtr(pd.wA);
   pd.wAt = pd.wA + szA;
   pd.wAt = ATL_AlignPtr(pd.wAt);

   pd.wC = pd.wAt + szA;
   pd.wC = ATL_AlignPtr(pd.wC);
   pd.A = A;
   pd.lda = lda;
   pd.C = C;
   pd.ldc = ldc;
   for (i=0; i < ndiag; i++)
   {
      pd.KbegCtr[i] = ATL_SetGlobalAtomicCount(nkcnt, nkblks, 0);
/*
 *    This blows up contention, but needed for correctness until we add
 *    ability to have last non-zero be 1!
 */
      pd.KdonCtr[i] = ATL_SetGlobalAtomicCountDown(nkcnt, nkblks);
      pd.Cdmuts[i] = ATL_mutex_init();
   }
   pd.cdmut = ATL_mutex_init();
   pd.cwmut = ATL_mutex_init();
   pd.cpydonBV = ATL_NewBV(ndiag);
   pd.dbetaBV = ATL_NewBV(ndiag);
   pd.cblkBV = ATL_NewBV(ncblks);
   pd.cbetaBV = ATL_NewBV(ncblks);

//   #define DEBUG1
   #ifdef DEBUG1
   {
      ATL_LAUNCHSTRUCT_t ls;
      ATL_thread_t ts;
      ts.rank = 0;
      ts.P = 1;
      ls.opstruct = &pd;
      Mjoin(PATL,DoWork_syrk_amm)(&ls, &ts);
   }
   #else
      ATL_goparallel(P, Mjoin(PATL,DoWork_syrk_amm), &pd, NULL);
   #endif
   #ifdef DEBUG2
      ATL_assert(ATL_FindFirstUnsetBitBV(pd.cblkBV, 0) == -1);
   #endif

/*
 * Free allocated structures and return;
 */
   ATL_FreeBV(pd.cpydonBV);
   ATL_FreeBV(pd.cblkBV);
   ATL_FreeBV(pd.dbetaBV);
   ATL_FreeBV(pd.cbetaBV);
   ATL_mutex_free(pd.cdmut);
   ATL_mutex_free(pd.cwmut);
   ATL_FreeGlobalAtomicCount(pd.dCtr);
   for (i=0; i < ndiag; i++)
   {
      ATL_FreeGlobalAtomicCount(pd.KbegCtr[i]);
      ATL_FreeGlobalAtomicCountDown(pd.KdonCtr[i]);
      ATL_mutex_free(pd.Cdmuts[i]);
   }
   free(vp);
   return(0);
}

void Mjoin(PATL,tsyrk_amm)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans, ATL_CINT N,
    ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const SCALAR beta, TYPE *C, ATL_CINT ldc)

{
   if (N <= Mmax(ATL_AMM_MAXMB,ATL_AMM_MAXNB))
   {
      if (!Mjoin(PATL,tsyrk_amm_K)(Uplo, Trans, N, K, alpha, A, lda,
                                   beta, C, ldc))
         return;
   }
/*
 *  Recur on K until tsyrk_amm can allocate enough space
 */
    if (Mjoin(PATL,tsyrk_amm_N)(Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc))
    {
       unsigned int kL = K>>1, kR=K-kL;
       const TYPE *a = (Trans == AtlasNoTrans) ? A+((size_t)lda)*kL : A+kL;
/*
 *     This stopping criteria should never happen, but it's here in case
 *     we have a system where you can't malloc much of anything, where we'll
 *     just try serial
 */
       if (kL < 32 || kR < 32)
       {
           Mjoin(PATL,syrk)(Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
           return;
       }
       Mjoin(PATL,tsyrk_amm)(Uplo, Trans, N, kL, alpha, A, lda, beta, C, ldc);
       Mjoin(PATL,tsyrk_amm)(Uplo, Trans, N, kR, alpha, a, lda,
                             ATL_rone, C, ldc);
    }
}
