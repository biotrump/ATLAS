#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_blk.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_ablk2cmat.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_cm2am_a1.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_cm2am_an.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_cm2am_aX.h))
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_flag.h))
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_kern.h))

static int IK=ATL_AMM_NCASES-1, MB=0, NB=0, KB=0;
#ifdef DCPLX
   static char MY_PRE='z', MY_PRE2='Z';
#elif defined(SCPLX)
   static char MY_PRE='c', MY_PRE2='C';
#elif defined(SREAL)
   static char MY_PRE='s', MY_PRE2='S';
#else
   static char MY_PRE='d', MY_PRE2='D';
#endif
/*
 * This routine overrides normal GetAmmmInfo, so we can tune kernel params
 */
int Mjoin(PATL,GetAmmmInfo)
(
   amminfo_t *out,
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alpha,
   const SCALAR beta
)
{
   int ik=IK, appAl;  /* 0:A, 1:B, 2:C */

   ATL_assert(ik >= 0 && ik < ATL_AMM_NCASES);
   while (K < ATL_AMM_KBs[ik] && ik)
      ik--;
   out->IDX = ik;
   out->mb = (MB) ? MB : ATL_AMM_MBs[ik];
   out->nb = (NB) ? NB : ATL_AMM_NBs[ik];
   out->kb = (KB) ? KB : ATL_AMM_KBs[ik];
   out->kbmin = ATL_AMM_KBMINs[ik];
   out->mu = ATL_AMM_MUs[ik];
   out->nu = ATL_AMM_NUs[ik];
   out->ku = ATL_AMM_KUs[ik];
   out->flag = ATL_AMM_KFLAG[ik];
   out->amm_b0 = ATL_AMM_KERN_b0[ik];
   out->amm_b1 = ATL_AMM_KERN_b1[ik];
   out->amm_bn = ATL_AMM_KERN_bn[ik];
   out->amm_k1_b0 = ATL_AMM_KERN_K1[ik];
   out->amm_k1_b1 = ATL_AMM_KERN_K1_b1[ik];
   out->amm_k1_bn = ATL_AMM_KERN_K1_bn[ik];
/*
 * Apply alpha to smallest matrix, and use alpha/beta to pick copy routines
 */
   if (SCALAR_IS_ONE(alpha))
   {
      appAl = 0;
      #ifdef TCPLX
         if (TA == AtlasNoTrans)
            out->a2blk = ATL_AMM_AT2BLK_a1[ik];
         else if (TA == AtlasTrans)
            out->a2blk = ATL_AMM_A2BLK_a1[ik];
         else if (TA == AtlasConjTrans)
            out->a2blk = ATL_AMM_AC2BLK_a1[ik];
         else
            out->a2blk = ATL_AMM_AH2BLK_a1[ik];
         if (TB == AtlasNoTrans)
             out->b2blk = ATL_AMM_B2BLK_a1[ik];
         else if (TB == AtlasTrans)
             out->b2blk = ATL_AMM_BT2BLK_a1[ik];
         else if (TB == AtlasConjTrans)
             out->b2blk = ATL_AMM_BH2BLK_a1[ik];
         else
             out->b2blk = ATL_AMM_BC2BLK_a1[ik];
      #else
         out->a2blk = (TA == AtlasNoTrans) ?
            ATL_AMM_AT2BLK_a1[ik]:ATL_AMM_A2BLK_a1[ik];
         out->b2blk = (TB == AtlasNoTrans) ?
            ATL_AMM_B2BLK_a1[ik]:ATL_AMM_BT2BLK_a1[ik];
      #endif
      if (SCALAR_IS_ONE(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_a1_b1[ik];
      else if (SCALAR_IS_ZERO(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_a1_b0[ik];
      else if (SCALAR_IS_NONE(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_a1_bn[ik];
      else
         out->Cblk2cm = ATL_AMM_BLK2C_a1_bX[ik];
   }
   else  /* alpha is not one */
   {
      if (M >= N)                  /* A is larger than B, put alpha on C or B */
         appAl = (M >= K) ? 1 : 2;
      else                         /* B is larger than A, put alpha on C or A */
         appAl = (N >= K) ? 0 : 2;
      if (appAl == 2)  /* apply alpha to C */
      {
         #ifdef TCPLX
            if (TA == AtlasNoTrans)
               out->a2blk = ATL_AMM_AT2BLK_a1[ik];
            else if (TA == AtlasTrans)
               out->a2blk = ATL_AMM_A2BLK_a1[ik];
            else if (TA == AtlasConjTrans)
               out->a2blk = ATL_AMM_AC2BLK_a1[ik];
            else
               out->a2blk = ATL_AMM_AH2BLK_a1[ik];
            if (TB == AtlasNoTrans)
                out->b2blk = ATL_AMM_B2BLK_a1[ik];
            else if (TB == AtlasTrans)
                out->b2blk = ATL_AMM_BT2BLK_a1[ik];
            else if (TB == AtlasConjTrans)
                out->b2blk = ATL_AMM_BH2BLK_a1[ik];
            else
                out->b2blk = ATL_AMM_BC2BLK_a1[ik];
         #else
            out->a2blk = (TA == AtlasNoTrans) ?
                         ATL_AMM_AT2BLK_a1[ik] : ATL_AMM_A2BLK_a1[ik];
            out->b2blk = (TB == AtlasNoTrans) ?
                         ATL_AMM_B2BLK_a1[ik] : ATL_AMM_BT2BLK_a1[ik];
         #endif
         if (SCALAR_IS_ONE(beta))
            out->Cblk2cm = SCALAR_IS_NONE(alpha) ?
                           ATL_AMM_BLK2C_an_b1[ik] : ATL_AMM_BLK2C_aX_b1[ik];
         else if (SCALAR_IS_ZERO(beta))
            out->Cblk2cm = SCALAR_IS_NONE(alpha) ?
                           ATL_AMM_BLK2C_an_b0[ik] : ATL_AMM_BLK2C_aX_b0[ik];
         else if (SCALAR_IS_NONE(beta))
            out->Cblk2cm = SCALAR_IS_NONE(alpha) ?
                           ATL_AMM_BLK2C_an_bn[ik] : ATL_AMM_BLK2C_aX_bn[ik];
         else
            out->Cblk2cm = SCALAR_IS_NONE(alpha) ?
                           ATL_AMM_BLK2C_an_bX[ik] : ATL_AMM_BLK2C_aX_bX[ik];
      }
      else  /* not applying alpha to C */
      {
         if (SCALAR_IS_ONE(beta))
            out->Cblk2cm = ATL_AMM_BLK2C_a1_b1[ik];
         else if (SCALAR_IS_ZERO(beta))
            out->Cblk2cm = ATL_AMM_BLK2C_a1_b0[ik];
         else if (SCALAR_IS_NONE(beta))
            out->Cblk2cm = ATL_AMM_BLK2C_a1_bn[ik];
         else
            out->Cblk2cm = ATL_AMM_BLK2C_a1_bX[ik];
         if (!appAl)  /* apply to alpha to A */
         {
            #ifdef TCPLX
               if (TB == AtlasNoTrans)
                   out->b2blk = ATL_AMM_B2BLK_a1[ik];
               else if (TB == AtlasTrans)
                   out->b2blk = ATL_AMM_BT2BLK_a1[ik];
               else if (TB == AtlasConjTrans)
                   out->b2blk = ATL_AMM_BH2BLK_a1[ik];
               else
                   out->b2blk = ATL_AMM_BC2BLK_a1[ik];
               if (SCALAR_IS_NONE(alpha))
               {
                  if (TA == AtlasNoTrans)
                     out->a2blk = ATL_AMM_AT2BLK_an[ik];
                  else if (TA == AtlasTrans)
                     out->a2blk = ATL_AMM_A2BLK_an[ik];
                  else if (TA == AtlasConjTrans)
                     out->a2blk = ATL_AMM_AC2BLK_an[ik];
                  else
                     out->a2blk = ATL_AMM_AH2BLK_an[ik];
               }
               else
               {
                  if (TA == AtlasNoTrans)
                     out->a2blk = ATL_AMM_AT2BLK_aX[ik];
                  else if (TA == AtlasTrans)
                     out->a2blk = ATL_AMM_A2BLK_aX[ik];
                  else if (TA == AtlasConjTrans)
                     out->a2blk = ATL_AMM_AC2BLK_aX[ik];
                  else
                     out->a2blk = ATL_AMM_AH2BLK_aX[ik];
               }
            #else
               if (SCALAR_IS_NONE(alpha))
                  out->a2blk = (TA == AtlasNoTrans) ?
                               ATL_AMM_AT2BLK_an[ik] : ATL_AMM_A2BLK_an[ik];
               else
                  out->a2blk = (TA == AtlasNoTrans) ?
                               ATL_AMM_AT2BLK_aX[ik] : ATL_AMM_A2BLK_aX[ik];
               out->b2blk = (TB == AtlasNoTrans) ?
                            ATL_AMM_B2BLK_a1[ik] : ATL_AMM_BT2BLK_a1[ik];
            #endif
         }
         else /* apply alpha to B */
         {
            #ifdef TCPLX
               if (TA == AtlasNoTrans)
                  out->a2blk = ATL_AMM_AT2BLK_a1[ik];
               else if (TA == AtlasTrans)
                  out->a2blk = ATL_AMM_A2BLK_a1[ik];
               else if (TA == AtlasConjTrans)
                  out->a2blk = ATL_AMM_AC2BLK_a1[ik];
               else
                  out->a2blk = ATL_AMM_AH2BLK_a1[ik];
               if (SCALAR_IS_NONE(alpha))
               {
                  if (TB == AtlasNoTrans)
                      out->b2blk = ATL_AMM_B2BLK_an[ik];
                  else if (TB == AtlasTrans)
                      out->b2blk = ATL_AMM_BT2BLK_an[ik];
                  else if (TB == AtlasConjTrans)
                      out->b2blk = ATL_AMM_BH2BLK_an[ik];
                  else
                      out->b2blk = ATL_AMM_BC2BLK_an[ik];
               }
               else
               {
                  if (TB == AtlasNoTrans)
                      out->b2blk = ATL_AMM_B2BLK_aX[ik];
                  else if (TB == AtlasTrans)
                      out->b2blk = ATL_AMM_BT2BLK_aX[ik];
                  else if (TB == AtlasConjTrans)
                      out->b2blk = ATL_AMM_BH2BLK_aX[ik];
                  else
                      out->b2blk = ATL_AMM_BC2BLK_aX[ik];
               }
            #else
               out->a2blk = (TA == AtlasNoTrans) ?
                            ATL_AMM_AT2BLK_a1[ik] : ATL_AMM_A2BLK_a1[ik];
               if (SCALAR_IS_NONE(alpha))
                  out->b2blk = (TB == AtlasNoTrans) ?
                               ATL_AMM_B2BLK_an[ik] : ATL_AMM_BT2BLK_an[ik];
               else
                  out->b2blk = (TB == AtlasNoTrans) ?
                               ATL_AMM_B2BLK_aX[ik] : ATL_AMM_BT2BLK_aX[ik];
            #endif
         }
      }
   }
   return(appAl);
}

double time00();
double Time2Mflops(int M, int N, int K, double t0)
{
   return((((2.0*M)*N)*K) / (1000000.0*t0));
}

int FindLowerBlock(int *M, int *N, int *K)
{
   int iret=1;
   if (ATL_AMM_KRUNTIME(ATL_AMM_KFLAG[IK]) && *K > *M && *K > *N &&
       *K-ATL_AMM_KUs[IK] >=  ATL_AMM_KBMINs[IK])
      *K -= ATL_AMM_KUs[IK];
   else if (*M >= *N && *M > ATL_AMM_MUs[IK])
      *M -= ATL_AMM_MUs[IK];
   else if (*N > ATL_AMM_NUs[IK])
      *N -= ATL_AMM_NUs[IK];
   else
      iret = 0;
   return(iret);
}

double GetMflops(ATL_CINT M, ATL_CINT N, ATL_CINT K, TYPE *A, TYPE *B, TYPE *C)
{
   double t0;
   const TYPE ONE[2] = {ATL_rone, ATL_rzero};

   t0 = time00();
   Mjoin(PATL,ammm)(AtlasNoTrans, AtlasNoTrans, M, N, K, ONE, A, M, B, K,
                    ONE, C, M);
   t0 = time00() - t0;
   return(Time2Mflops(M, N, K, t0));
}

void TuneBlocking(int M, int N, int K, TYPE *A, TYPE *B, TYPE *C,
              int *MB_, int *NB_, int *KB_)
{
   int mbB=ATL_AMM_MBs[IK], nbB=ATL_AMM_NBs[IK], kbB=ATL_AMM_KBs[IK];
   int mb=mbB, nb=nbB, kb=kbB, m, n, k;
   double t0, t1, mf0, mf, mfB;

   printf("\nTRYING REDUCED BLOCKING WITH INDEX=%d\n", IK);
   printf("        M       N       K  IDX   MB   NB   KB         MFLOPS\n");
   printf("   ======  ======  ======  ===  ===  ===  ===  =============\n");
   m = (M/mbB)*mbB;
   n = (N/nbB)*nbB;
   k = (K/kbB)*kbB;
   mf0 = mf = mfB = GetMflops(m, n, k, A, B, C);
   printf("  %7d %7d %7d %4d %4d %4d %4d %14.2f\n",
          m, n, k, IK, mb, nb, kb, mf);
   while(FindLowerBlock(&mb, &nb, &kb))
   {
      MB = mb;
      NB = nb;
      KB = kb;
      m = (M/mb)*mb;
      n = (N/nb)*nb;
      k = (K/kb)*kb;
      mf = GetMflops(m, n, k, A, B, C);
      printf("  %7d %7d %7d %4d %4d %4d %4d %14.2f\n",
             m, n, k, IK, mb, nb, kb, mf);
      if (mf > mfB)
      {
         mbB = mb;
         nbB = nb;
         kbB = kb;
         mfB = mf;
      }
      else if (mf*1.02 <= mfB)
         break;
   }
   printf("BEST BLOCKING: MB=%d, NB=%d, KB=%d speedup=%.3f\n",
          mbB, nbB, kbB, (mfB/mf0));
   *MB_ = mbB;
   *NB_ = nbB;
   *KB_ = kbB;
}

int TuneIndx(int M, int N, int K, TYPE *A, TYPE *B, TYPE *C)
{
   double mfB=0.0, mf;
   int ik, iB=ATL_AMM_NCASES-1;

   printf("\n\nFINDING BEST INDEX:\n");
   printf("        M       N       K  IDX   MB   NB   KB         MFLOPS\n");
   printf("   ======  ======  ======  ===  ===  ===  ===  =============\n");
   for (ik=iB; ik >= 0; ik--)
   {
      const int mb=ATL_AMM_MBs[ik],nb=ATL_AMM_NBs[ik],kb=ATL_AMM_KBs[ik];
      const int m=(M/mb)*mb, n=(N/nb)*nb, k=(K/kb)*kb;
      IK = ik;
      mf = GetMflops(m, n, k, A, B, C);
      printf("  %7d %7d %7d %4d %4d %4d %4d %14.2f\n",
             m, n, k, ik, mb, nb, kb, mf);
      if (mf > mfB)
      {
         iB = ik;
         mfB = mf;
      }
      else if (ATL_AMM_KBs[ik] < 16)
         break;
   }
   printf("\nBEST INDEX=%d (%.2f)\n", iB, mfB);
   return(iB);
}

void TuneCplxNB(FILE *fpout, int M, int N, int K, TYPE *A, TYPE *B, TYPE *C)
{
   int ik, mb, nb, kb;
   ik = TuneIndx(M, N, K, A, B, C);
   IK = ik;
   TuneBlocking(M, N, K, A, B, C, &mb, &nb, &kb);
   fprintf(fpout, "#ifndef ATL_%cAMM_SUM_H\n   #define ATLAS_%cAMM_SUM_H\n\n",
           MY_PRE2, MY_PRE2);
   fprintf(fpout, "   #define ATL_CAMM_MAXINDX %d\n", IK);
   fprintf(fpout, "   #define ATL_CAMM_MAXMB %d\n", mb);
   fprintf(fpout, "   #define ATL_CAMM_MAXNB %d\n", nb);
   fprintf(fpout, "   #define ATL_CAMM_MAXKB %d\n", kb);
   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
}

FILE *GetFlags(int nargs, char **args, int *M, int *N, int *K)
{
   FILE *fp=NULL;
   *K = 1200;
   *M = *N = 2000;
   if (!fp)
   {
      char nam[32];
      sprintf(nam, "res/atlas_%samm_sum.h", Mstr(PRE));
      fp = fopen(nam, "w");
      ATL_assert(fp);
   }
   return(fp);
}


int main(int nargs, char **args)
{
   void *vp;
   TYPE *A, *B, *C;
   int M, N, K;
   size_t szA, szB, szC;
   FILE *fpout;

   fpout = GetFlags(nargs, args, &M, &N, &K);
   M = ((M+ATL_AMM_MAXMB-1)/ATL_AMM_MAXMB)*ATL_AMM_MAXMB;
   N = ((N+ATL_AMM_MAXNB-1)/ATL_AMM_MAXNB)*ATL_AMM_MAXNB;
   K = ((K+ATL_AMM_MAXKB-1)/ATL_AMM_MAXKB)*ATL_AMM_MAXKB;
   szA = M*K;
   szB = K*N;
   szC = M*N;
   vp = malloc(ATL_MulBySize(szA + szB + szC) + 3*ATL_Cachelen);
   ATL_assert(vp);
   A = ATL_AlignPtr(vp);
   B = A + szA+szA;
   B = ATL_AlignPtr(B);
   C = B + szB+szB;
   C = ATL_AlignPtr(C);
   Mjoin(PATL,zero)(szA, A, 1);
   Mjoin(PATL,zero)(szB, B, 1);
   Mjoin(PATL,zero)(szC, C, 1);
   TuneCplxNB(fpout, M, N, K, A, B, C);
   free(vp);
   return(0);
}
