#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"

static void DoCpy
(
   ATL_tgemm_ammG_t *pd,             /* problem definition */
   unsigned const int rank,          /* my (virtual) rank */
   unsigned const int dblk,          /* CpyBlk of C to compute */
   const size_t i,                   /* row index to start at */
   const size_t j,                   /* col index to start at */
   unsigned int kblk,                /* kblk to compute */
   unsigned const int mb,            /* actual sz of M w/o padding */
   unsigned const int nb             /* actual sz of N w/o padding */
)
{
   const int CPYA = dblk < pd->nmblks, CPYB = dblk < pd->nnblks;
   int kb = pd->kb, KB0=kb;
   const TYPE *A;
   const TYPE *B;
   TYPE *wA, *wB;
   size_t k;

   if (kblk)
      k = pd->kb0 + (kblk-1)*kb;
   else
   {
      k = 0;
      if (kb != pd->kb0)
      {
         kb = pd->kb0;
         KB0 = pd->KB0;
      }
   }

   if (CPYA)
   {
      A = pd->A + ((pd->TA) ? k + i*pd->lda : i + k*pd->lda);
      wA = pd->wA + kblk*pd->blkszA;
      wA += pd->panszA * ((CPYA) ? dblk : (pd->nmblks-1));
      pd->a2blk(kb, mb, pd->alpA, A, pd->lda, wA);
   }
   if (CPYB)
   {
      B = pd->B + ((pd->TB) ?  j + k*pd->ldb : k + j*pd->ldb);
      wB = pd->wB + kblk*pd->blkszB;
      wB += pd->panszB * ((CPYB) ? dblk : (pd->nnblks-1));
      pd->b2blk(kb, nb, pd->alpB, B, pd->ldb, wB);
   }
/*
 * If I just finished copying the K-panel for A or B, let folks know
 */
   if (ATL_DecGlobalAtomicCountDown((pd->KdonCtr)[dblk], rank) == 1)
   {
      ATL_mutex_lock(pd->cpmut);
      if (CPYA)
      {
         ATL_SetBitBV(pd->cpyAdBV, dblk);
         if (ATL_FindFirstUnsetBitBV(pd->cpyAdBV, 0) == -1)
            pd->cpyAdone = 1;
      }
      if (CPYB)
      {
         ATL_SetBitBV(pd->cpyBdBV, dblk);
         if (ATL_FindFirstUnsetBitBV(pd->cpyBdBV, 0) == -1)
            pd->cpyBdone = 1;
      }
      ATL_mutex_unlock(pd->cpmut);
   }
}

static void DoBlkWtCpy
(
   ATL_tgemm_ammG_t *pd,             /* problem definition */
   unsigned const int rank,          /* my (virtual) rank */
   unsigned const int dblk,          /* CpyBlk of C to compute */
   const size_t i,                   /* row index to start at */
   const size_t j,                   /* col index to start at */
   unsigned int kblk,                /* kblk to compute */
   unsigned const int mb,            /* actual sz of M w/o padding */
   unsigned const int nb,            /* actual sz of N w/o padding */
   unsigned const int nmu,           /* CEIL(mb/mu) */
   unsigned const int nnu,           /* CEIL(nb/nu) */
   ammkern_t amm,                    /* kern to use for all non-kb0 blocks */
   TYPE *wC                          /* my private C workspace */
)
{
   int kb = pd->kb, KB0=kb;
   const TYPE *A = pd->A;
   const TYPE *B = pd->B;
   TYPE *wA, *wB;
   size_t k;

   if (kblk)
      k = pd->kb0 + (kblk-1)*kb;
   else
   {
      k = 0;
      if (kb != pd->kb0)
      {
         kb = pd->kb0;
         KB0 = pd->KB0;
         if (amm == pd->amm_b0)
            amm = pd->ammK_b0;
         else
            amm = pd->ammK_b1;
      }
   }

   A += (pd->TA) ? k + i*pd->lda : i + k*pd->lda;
   B += (pd->TB) ?  j + k*pd->ldb : k + j*pd->ldb;
   wA = pd->wA + kblk*pd->blkszA + pd->panszA * dblk;
   wB = pd->wB + kblk*pd->blkszB + pd->panszB * dblk;
   pd->a2blk(kb, mb, pd->alpA, A, pd->lda, wA);
   pd->b2blk(kb, nb, pd->alpB, B, pd->ldb, wB);
/*
 * If I just finished copying the K-panel for both A & B, let folks know
 */
   if (ATL_DecGlobalAtomicCountDown((pd->KdonCtr)[dblk], rank) == 1)
   {
      ATL_mutex_lock(pd->cpmut);
      ATL_SetBitBV(pd->cpyAdBV, dblk);
      if (ATL_FindFirstUnsetBitBV(pd->cpyAdBV, 0) == -1)
         pd->cpyAdone = 1;

      ATL_SetBitBV(pd->cpyBdBV, dblk);
      if (ATL_FindFirstUnsetBitBV(pd->cpyBdBV, 0) == -1)
         pd->cpyBdone = 1;
      ATL_mutex_unlock(pd->cpmut);
   }
   amm(nmu, nnu, KB0, wA, wB, wC, wA, wB, wC);
}
/*
 * This routine computes one block of C, and writes answer back out to C
 */
static void DoBlksWtCpy
(
   ATL_tgemm_ammG_t *pd,             /* problem definition */
   unsigned const int rank,          /* my (virtual) rank */
   unsigned const int dblk,          /* CpyBlk of C to compute */
   unsigned int kctr,                /* non-zero KbegCtr */
   TYPE *wC                          /* my private C workspace */
)
{
   const int nmblks=pd->nmblks, nnblks=pd->nnblks, nkblks=pd->nkblks;
   const int minblks=Mmin(nmblks, nnblks);
   const int mb = (dblk >= nmblks-1) ? pd->mbf : pd->mb;
   const int nb = (dblk >= nnblks-1) ? pd->nbf : pd->nb;
   const int nmu = (dblk >= nmblks-1) ? pd->nmuf : pd->nmu;
   const int nnu = (dblk >= nnblks-1) ? pd->nnuf : pd->nnu;
   size_t i, j;
   unsigned int kblk = pd->nkblks - kctr;
   TYPE *wA, *wB, *c;
   TYPE beta = ATL_rone;

   i = (dblk < nmblks) ? dblk : nmblks-1;
   j = (dblk < nnblks) ? dblk : nnblks-1;
//fprintf(stderr, "%d: (%d,%d)\n", __LINE__, i, j);
   i *= pd->mb;
   j *= pd->nb;
   if (dblk < minblks)
      DoBlkWtCpy(pd, rank, dblk, i, j, kblk, mb, nb, nmu, nnu, pd->amm_b0, wC);
   else
      DoCpy(pd, rank, dblk, i, j, kblk, mb, nb);
   while (kctr = ATL_DecGlobalAtomicCount((pd->KbegCtr)[dblk], rank))
   {
      kblk = nkblks - kctr;
      if (dblk < minblks)
         DoBlkWtCpy(pd, rank, dblk, i, j, kblk, mb, nb, nmu, nnu,
                    pd->amm_b1, wC);
      else
         DoCpy(pd, rank, dblk, i, j, kblk, mb, nb);
   }
/*
 * Seize mutex for block of original C, and add my part out
 */
   if (dblk < minblks)
   {
      int COPIED=0;
      c = pd->C + j*pd->ldc + i;
      ATL_mutex_lock(pd->Cmuts[dblk]);
      if (!ATL_IsBitSetBV(pd->cbetaBV, dblk))
      {
         ATL_mutex_lock(pd->cbetamut);
         if (!ATL_IsBitSetBV(pd->cbetaBV, dblk))
         {
            ATL_SetBitBV(pd->cbetaBV, dblk);
            ATL_mutex_unlock(pd->cbetamut);
            pd->blk2c(mb, nb, pd->alpC, wC, pd->beta, c, pd->ldc);
            COPIED = 1;
         }
         else
            ATL_mutex_unlock(pd->cbetamut);
      }
      if (!COPIED)
         pd->blk2c_b1(mb, nb, pd->alpC, wC, ATL_rone, c, pd->ldc);
      ATL_mutex_unlock(pd->Cmuts[dblk]);
   }
}

static void DoCopyBlksG(ATL_tgemm_ammG_t *pd, int rank, int vrank, TYPE *wC)
{
   int NEWBLK=1;
   const int ND = pd->nMNblks;
/*
 * Keep going as long as there is copy work to be done
 */
   while (!pd->NOCPWORK)
   {
       int d=0, k;
/*
 *     Find which diagonal block to work on, and then which k blk to use
 */
       if (NEWBLK)
       {
          d = ATL_DecGlobalAtomicCount(pd->ccCtr, vrank);
          if (d)
          {
             k = ATL_DecGlobalAtomicCount((pd->KbegCtr)[ND-d], rank);
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
          unsigned int i;
          NEWBLK = 0;
          for (i=0; i < ND; i++)
          {
             unsigned int j = (i+vrank)%ND;
             k = ATL_DecGlobalAtomicCount((pd->KbegCtr)[j], vrank);
             d = ND-j;
             if (k)
                goto FOUNDDK;
          }
          pd->NOCPWORK = 1;   /* no work left to assign */
          return;             /* so done with copy of A&B and this func */
       }
/*
 *     If I reach here, I've got a valid d & k;  and I'll call a routine
 *     that continues to grab blocks from this diag along K until all K
 *     is done; it will then write the answer back to the original C, and
 *     return to this loop to see if it can help with another diag.
 */
       FOUNDDK:
          DoBlksWtCpy(pd, vrank, ND-d, k, wC);
   }
}

/*
 * RETURNS: undone C block that has had its A and B k-panels copied or:
 *          -1: no work left to do;     -2: no work available
 * NOTE: this routine only called when A/B copying is still ongoing!
 */
static int FindCblk(const int rank, ATL_tgemm_ammG_t *pd)
{
   unsigned const int ncblks=pd->nCblks, nmblks=pd->nmblks, nnblks=pd->nnblks;
   int bb=0, ib;
   do
   {
      int i, j;
      ib = ATL_FindFirstUnsetBitBV(pd->cCblkBV, bb);
      if (ib == -1)
         return(-1);
      j = ib / nmblks;
      i = ib - j*nmblks;
      if (pd->cpyAdone || ATL_IsBitSetBV(pd->cpyAdBV, i))
      {
         if (pd->cpyBdone || ATL_IsBitSetBV(pd->cpyBdBV, j))
            return(ib);
      }
      bb = ib+1;
   }
   while(bb < ncblks);
   return(-2);
}

/*
 * RETURNS: C block that has had its A and B k-panels copied
 */
static INLINE int GetCblk(const int rank, ATL_tgemm_ammG_t *pd)
{
   const int ncblks = pd->nCblks;
   int ib;
/*
 * If we are done copying, then use the counter to get the C block to work on
 * NOTE: thread that sets cpy[A,B]done must have cpmut to avoid race cond!
 */
   if (pd->cpyAdone & pd->cpyBdone)
   {
      do
      {
         ib = ATL_DecGlobalAtomicCount(pd->cCtr, rank);
         if (!ib)
            return(-1);
         ib = ncblks - ib;
      }
      while (ATL_IsBitSetBV(pd->cCblkBV, ib));
      return(ib);
   }
/*
 * If we reach here, we must assign C work based on dep info in cpy[A,B]dBV
 */
   else
   {
      while(1)
      {
         ib = FindCblk(rank, pd);
/*
 *       If we have a candidate block, grab the mutex and make sure one is
 *       still there
 */
         if (ib >= 0)
         {
            ATL_mutex_lock(pd->cpmut);
/*
 *          If copy finished while locking, call ourselves to get fast answer
 */
            if (pd->cpyAdone & pd->cpyBdone)
            {
               ATL_mutex_unlock(pd->cpmut);
               return(GetCblk(rank, pd));
            }
/*
 *          Otherwise, see if we've still got a candidate block
 */
            ib = FindCblk(rank, pd);
            if (ib >= 0)  /* we've got a block! */
            {
               ATL_SetBitBV(pd->cCblkBV, ib);  /* claim the block */
               ATL_mutex_unlock(pd->cpmut);
               return(ib);
            }
            ATL_mutex_unlock(pd->cpmut);
         }
         if (ib == -1)
            return(-1);
/*
 *       If no work available, pause before trying again
 */
         if (ib == -2)
            ATL_thread_yield();
      }
   }
}

/*
 * Computes (i,j) block of C whose data has previously been copied
 */
static void DoCblk(ATL_tgemm_ammG_t *pd, int rank, TYPE *wC, int i, int j)
{
   const unsigned int szA=pd->blkszA, szB=pd->blkszB, nkblks=pd->nkblks,
      kb=pd->kb;
   unsigned int nmu, nnu, mb, nb, k;
   const TYPE *wA = pd->wA + i*pd->panszA;
   const TYPE *wB = pd->wB + j*pd->panszB;
   const TYPE *wAn = wA+szA, *wBn = wB+szB;
   const ammkern_t amm = pd->amm_b1;

   #ifdef DEBUG2
     /* fprintf(stderr, "%d: DoCblk, (%d,%d)\n", __LINE__, i, j); */
      ATL_assert(ATL_IsBitSetBV(pd->cpyAdBV, i));
      ATL_assert(ATL_IsBitSetBV(pd->cpyBdBV, j));
   #endif
   if (i != pd->nmblks-1)
   {
      nmu = pd->nmu;
      mb = pd->mb;
   }
   else
   {
      nmu = pd->nmuf;
      mb = pd->mbf;
   }
   if (j != pd->nnblks-1)
   {
      nnu = pd->nnu;
      nb = pd->nb;
   }
   else
   {
      nnu = pd->nnuf;
      nb = pd->nbf;
   }
   pd->ammK_b0(nmu, nnu, pd->KB0, wA, wB, wC, wAn, wBn, wC);
   for (k=1; k < nkblks; k++)
   {
      wA = wAn;
      wB = wBn;
      wAn += szA;
      wBn += szB;
      amm(nmu, nnu, kb, wA, wB, wC, wAn, wBn, wC);
   }
   pd->blk2c(mb, nb, pd->alpC, wC, pd->beta, pd->C+pd->ldc*j*pd->nb+i*pd->mb,
             pd->ldc);
}

static void DoNoCopyBlksG(ATL_tgemm_ammG_t *pd, int rank, int vrank, TYPE *wC)
{
   const unsigned int nmblks = pd->nmblks;
   int ic;
   while ((ic = GetCblk(rank, pd)) != -1)
   {
      int i, j;
      j = ic / nmblks;
      i = ic - j*nmblks;
      DoCblk(pd, vrank, wC, i, j);
   }
}

void Mjoin(PATL,DoWork_ammG)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tgemm_ammG_t *pd = pp->PD;
   TYPE *myC = pd->wC + vrank*pd->blkszC;
   if (!pd->NOCPWORK)
      DoCopyBlksG(pd, rank, vrank, myC);
   DoNoCopyBlksG(pd, rank, vrank, myC);
}
/*
 * This one can be used anytime C is large enough to provide parallelism.
 * It tries to copy all of A & B up front, so recursion may be needed to
 * use it for large matrices.
 * RETURNS: 0 if it did the operation, non-zero if it did not
 */
int Mjoin(PATL,tammm_G)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CINT M,
   ATL_CINT N,
   ATL_CINT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CINT lda,
   const TYPE *B,
   ATL_CINT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CINT ldc
)
{
   int ncblks, nmblks, nnblks, nkblks, nMNblks, mb, nb, kb, mu, nu, ku;
   int kb0, KB0, P, nkcnt, i, j, k;
   size_t sz, szA, szB;
   ATL_tgemm_ammG_t pd;
   amminfo_t mminfo;
   void *vp;
   #ifdef ATL_PHI_SORTED
      const unsigned int nthr = ATL_TP_PTR ? ATL_TP_PTR->nthr : ATL_NTHREADS,
                         p4 = nthr>>2, p4_2 = nthr+nthr, p4_3 = p4_2+p4;
   #endif

   nmblks = M / ATL_AMM_66MB;
   nnblks = N / ATL_AMM_66NB;
   ncblks = sz = ((size_t)nmblks) * nnblks;
   if (ncblks != sz)   /* too many blocks to safely count */
      return(1);       /* so tell caller to recur or call something else */
/*
 * Quick exit for problems too small to thread
 */
   if (ncblks < 8)
      return(Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                              beta, C, ldc));
   mb = Mjoin(PATL,tGetAmmmInfo)(&mminfo, Mmin(ncblks, ATL_NTHREADS), TA,
                                 TB, M, N, K, alpha, beta);
   pd.alpA = pd.alpB = pd.alpC = ATL_rone;
   if (!mb)
      pd.alpA = alpha;
   else if (mb == 1)
      pd.alpB = alpha;
   else
      pd.alpC = alpha;
   pd.beta = beta;
   pd.blk2c = mminfo.Cblk2cm;
   pd.blk2c_b1 = mminfo.Cblk2cm_b1;
   pd.a2blk = mminfo.a2blk;
   pd.b2blk = mminfo.b2blk;
   pd.mb = mb = mminfo.mb;
   pd.nb = nb = mminfo.nb;
   pd.kb = kb = mminfo.kb;
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   pd.nmu = mb/mu;
   pd.nnu = nb/nu;
   pd.nmblks = nmblks = (M+mb-1)/mb;
   pd.nnblks = nnblks = (N+nb-1)/nb;
   pd.nkblks = nkblks = (K+kb-1)/kb;
   pd.nCblks = sz = ((size_t)nmblks)*nnblks;
   if (pd.nCblks != sz)
      return(1);
   P = Mmin(ATL_NTHREADS,pd.nCblks);
   pd.nMNblks = nMNblks = Mmax(nmblks, nnblks);
   sz = nmblks * mb;
   pd.mbf = (sz == M) ? mb : M+mb-sz;
   sz = nnblks * nb;
   pd.nbf = (sz == N) ? nb : N+nb-sz;
   pd.nmuf = (pd.mbf+mu-1)/mu;
   pd.nnuf = (pd.nbf+nu-1)/nu;
   pd.blkszA = mb * kb;
   pd.blkszB = nb * kb;
   pd.blkszC = mb * nb;
   KB0 = kb0 = K - (K/kb)*kb;
   if (!kb0)
   {
      kb0 = KB0 = kb;
      pd.ammK_b0 = mminfo.amm_b0;
      pd.ammK_b1 = mminfo.amm_b1;
   }
   else
   {
      #if ATL_AMM_MAXKMAJ > 1
         if (ATL_AMMFLG_KMAJOR(mminfo.flag))
         {
            KB0 = ((kb0+ku-1)/ku)*ku;
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
         if (ATL_AMMFLG_KRUNTIME(mminfo.flag) && kb0 == (kb0/ku)*ku &&
             kb0 > mminfo.kbmin)
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
   pd.kb0 = kb0;
   pd.KB0 = KB0;
   pd.amm_b0 = mminfo.amm_b0;
   pd.amm_b1 = mminfo.amm_b1;
   szA = nmblks*nkblks*pd.blkszA;
   szB = nnblks*nkblks*pd.blkszB;
   sz = nMNblks*3*sizeof(void*);    /* for K[beg,don]Ctr & Cmuts arrays */
   sz += ATL_MulBySize(pd.blkszC)*P;   /* local C wrkspc */
   sz += ATL_MulBySize(szA+szB);              /* A/B workspace */
   sz += 3*ATL_Cachelen;                      /* room for alignment */
   #ifdef ATL_PHI_SORTED
      pd.ncntxts = 1;
      pd.ncores = p4;
      if (P == p4_2)
         pd.ncntxts = 2;
      else if (P == p4_3)
         pd.ncntxts = 3;
      else if (P == nthr)
         pd.ncntxts = 4;
      else
         pd.ncores = P;
      if (pd.ncntxts > 1)
         sz += ATL_Cachelen*(pd.ncores+1);
   #endif
   if ((sz>>20) > ATL_PTMAXMALLOC_MB)
      return(2);
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
         *ip = ip[1] = ip[2] = ip[3] = -2;
      }
   }
   else
   #endif
      pd.KbegCtr = vp;
   pd.cpyAdone = pd.cpyBdone = pd.NOCPWORK = 0;
   pd.TA = (TA == AtlasTrans);
   pd.TB = (TB == AtlasTrans);
   pd.panszA = nkblks*pd.blkszA;
   pd.panszB = nkblks*pd.blkszB;
   pd.ccCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nMNblks, P),nMNblks, 0);
   pd.cCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nCblks, P),pd.nCblks, 0);
   nkcnt = ATL_EstNctr(nkblks, P);
   pd.KdonCtr = pd.KbegCtr + nMNblks;
   pd.Cmuts = pd.KdonCtr + nMNblks;
   pd.wA = (TYPE*)(pd.Cmuts+nMNblks);
   pd.wA = ATL_AlignPtr(pd.wA);
   pd.wB = pd.wA + szA;
   pd.wB = ATL_AlignPtr(pd.wB);
   pd.wC = pd.wB + szB;
   pd.wC = ATL_AlignPtr(pd.wC);
   pd.A = A;
   pd.B = B;
   pd.lda = lda;
   pd.ldb = ldb;
   pd.C = C;
   pd.ldc = ldc;
   for (i=0; i < nMNblks; i++)
   {
      pd.KbegCtr[i] = ATL_SetGlobalAtomicCount(nkcnt, nkblks, 0);
/*
 *    This blows up contention, but needed for correctness until we add
 *    ability to have last non-zero be 1!
 */
      pd.KdonCtr[i] = ATL_SetGlobalAtomicCountDown(nkcnt, nkblks);
      pd.Cmuts[i] = ATL_mutex_init();
   }
   pd.cpmut = ATL_mutex_init();
   pd.cbetamut = ATL_mutex_init();
   pd.cpyAdBV = ATL_NewBV(nmblks);
   pd.cpyBdBV = ATL_NewBV(nnblks);
   pd.cCblkBV = ATL_NewBV(pd.nCblks);
   pd.cbetaBV = ATL_NewBV(nMNblks);
/*
 * Initialize cCblkBV so that all diagonal blocks are shown as already complete
 * This BV is used for assigning work to non-copy blocks.
 */
   k = Mmin(nmblks, nnblks);
   for (i=0; i < k; i++)
      ATL_SetBitBV(pd.cCblkBV, i*nmblks+i);
   k = Mmin(nmblks, nnblks);

//   #define DEBUG1
   #ifdef DEBUG1
   {
      ATL_tpool_t *pp=ATL_TP_PTR;
      if (!pp)
         pp = ATL_NewThreadPool(1, 0, NULL);
      ATL_assert(pp);
      pp->PD = &pd;
      Mjoin(PATL,DoWork_ammG)(pp, 0, 0);
      if (pp != ATL_TP_PTR)
         ATL_FreeThreadPool(pp);
   }
   #else
      ATL_goParallel(P, Mjoin(PATL,DoWork_ammG), NULL, &pd, NULL);
   #endif
/*
 * Free allocated structures and return;
 */
   ATL_FreeBV(pd.cpyAdBV);
   ATL_FreeBV(pd.cpyBdBV);
   ATL_FreeBV(pd.cCblkBV);
   ATL_FreeBV(pd.cbetaBV);
   ATL_mutex_free(pd.cpmut);
   ATL_mutex_free(pd.cbetamut);
   ATL_FreeGlobalAtomicCount(pd.cCtr);
   ATL_FreeGlobalAtomicCount(pd.ccCtr);
   for (i=0; i < nMNblks; i++)
   {
      ATL_FreeGlobalAtomicCount(pd.KbegCtr[i]);
      ATL_FreeGlobalAtomicCountDown(pd.KdonCtr[i]);
      ATL_mutex_free(pd.Cmuts[i]);
   }
   free(vp);
   return(0);
}
