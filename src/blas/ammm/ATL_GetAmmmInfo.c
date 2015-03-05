#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_blk.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_ablk2cmat.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_cm2am_a1.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_cm2am_an.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_cm2am_aX.h))
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_flag.h))
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_kern.h))

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
   #ifdef ATL_CAMM_MAXINDX
      int ik=ATL_CAMM_MAXINDX;
   #else
      int ik=ATL_AMM_NCASES-1;
   #endif
   int appAl;  /* 0:A, 1:B, 2:C */
/*
 * For rank-K update, choose smallest KB that contains required K
 */
   if (K < ATL_AMM_MAXKB)
   {
      for (ik=0; ik < ATL_AMM_NCASES && ATL_AMM_KBs[ik] < K; ik++);
   }
   out->IDX = ik;
   out->mb = ATL_AMM_MBs[ik];
   out->nb = ATL_AMM_NBs[ik];
   out->kb = ATL_AMM_KBs[ik];
   #ifdef ATL_CAMM_MAXINDX
      if (ik == ATL_CAMM_MAXINDX)
      {
         out->mb = ATL_CAMM_MAXMB;
         out->nb = ATL_CAMM_MAXNB;
         out->kb = ATL_CAMM_MAXKB;
      }
   #endif
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

/*
 * Following routines help pick NB for parallel routines
 */
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_sum.h))
#ifdef ATL_CAMM_MAXINDX
   #define ATL_MAXIDX ATL_CAMM_MAXINDX
#elif defined(ATL_AMM_98IDX) && !defined(ATL_AMM_66IDX)
   #define ATL_MAXIDX ATL_AMM_98IDX
#endif
#ifndef ATL_MAXIDX
   #define ATL_MAXIDX ATL_AMM_NCASES-1
#endif


int Mjoin(PATL,tGetAmmmInfo)
(
   amminfo_t *out,
   const unsigned int P,
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alpha,
   const SCALAR beta
)
{
   int ik=ATL_MAXIDX, mb, nb, nmblks, nnblks, ncblks;
   int appAl;  /* 0:A, 1:B, 2:C */
/*
 * For rank-K update, choose smallest KB that contains required K
 */
   if (K < ATL_AMM_MAXKB)
   {
      for (ik=0; ik < ATL_AMM_NCASES && ATL_AMM_KBs[ik] < K; ik++);
   }
   ik++;
   do
   {
      ik--;
      mb = ATL_AMM_MBs[ik];
      nb = ATL_AMM_NBs[ik];
      nmblks = M / mb;
      nnblks = N / nb;
      ncblks = nmblks * nnblks;
   }
   while (ik >= ATL_AMM_66IDX && ncblks < P);
   out->IDX = ik;
   out->mb = mb;
   out->nb = nb;
   out->kb = ATL_AMM_KBs[ik];
   #ifdef ATL_CAMM_MAXINDX
      if (ik == ATL_CAMM_MAXINDX)
      {
         out->mb = ATL_CAMM_MAXMB;
         out->nb = ATL_CAMM_MAXNB;
         out->kb = ATL_CAMM_MAXKB;
      }
   #endif
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
      out->Cblk2cm_b1 = ATL_AMM_BLK2C_a1_b1[ik];
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
         out->Cblk2cm_b1 = SCALAR_IS_NONE(alpha) ?
                        ATL_AMM_BLK2C_an_b1[ik] : ATL_AMM_BLK2C_aX_b1[ik];
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
         out->Cblk2cm_b1 = ATL_AMM_BLK2C_a1_b1[ik];
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
#ifdef TREAL
/*
 * For TRSM, we're going to get our main parellelism from the N dimension
 * We know that mb==kb, but nb is independent.  alpha for A will be -1,
 * and for B it will 1.  Cblk2cm will be beta=0, while Cbk2cm_b1 will be 1.
 * RETURNS: upper bound on useful nthreads to use
 */
#include "atlas_ttypes.h"
int Mjoin(PATL,tGetTrsmInfo)
(
   ATL_ttrsm_amm_t *pd,
   int P,
   enum ATLAS_TRANS TA,
   ATL_CSZT M,
   ATL_CSZT N,
   const SCALAR beta
)
{
   #ifdef ATL_CAMM_MAXINDX
      int ik = ATL_CAMM_MAXINDX;
   #else
      int ik = ATL_AMM_NCASES-1;
   #endif
   int mu, nu, ku, nb, nnblks, mb, nmblks, mb0;

/*
 * Get a KB smaller than M
 */
   for (; ik > 0 && ATL_AMM_KBs[ik] > M; ik--);

/*
 * Find a kernel where mb can be set to kb; we know these exist, since we insist
 * on square problems of moderate size
 */
   for (; ik > 0; ik--)
   {
      mu = ATL_AMM_MUs[ik];
      mb = ATL_AMM_KBs[ik];
      if (mb > M)
         continue;
/*
 *    Any kernel can be used if it can be called with MB = KB
 */
      if ((mb/mu)*mu == mb)       /* it is legal to call wt MB=KB */
         break;                   /* so use this kernel */
/*
 *    KRUNTIME kernels can vary their KB, and thus be made legal
 */
      if (ATL_AMM_KRUNTIME(ATL_AMM_KFLAG[ik]))
      {
         ku = ATL_AMM_KUs[ik];
         ku = ATL_lcm(ku, mu);
         mb = (mb/ku)*ku;
         if (mb)
            break;
      }
   }
// ik=0;   /* FOR TESTING!!!!! */
   pd->mb = mb = ATL_AMM_KBs[ik];
   nb = ATL_AMM_NBs[ik];
   nu = ATL_AMM_NUs[ik];
   if (P*nb > N)
      P = (N+nb-1)/nb;
   pd->nb = nb;
   mu = ATL_AMM_MUs[ik];
   ku = ATL_AMM_KUs[ik];
   pd->nmu = mb / mu;
   pd->nnu = nb / nu;
   nnblks = N / nb;
   pd->nbf = N - nb*nnblks;
   if (!pd->nbf)
   {
      pd->nbf = nb;
      pd->nnuf = pd->nnu;
   }
   else
   {
      nnblks++;
      pd->nnuf = (pd->nbf+nu-1)/nu;
   }
   pd->nnblks = nnblks;
   pd->amm_b0 = ATL_AMM_KERN_b0[ik];
   pd->amm_b1 = ATL_AMM_KERN_b1[ik];
   nmblks = M/mb;
   mb0 = (M - nmblks*mb);
   if (!mb0)
   {
      pd->MB0 = mb0 = mb;
      pd->nmu0 = pd->nmu;
   }
   else
   {
      nmblks++;
      if (ATL_AMM_KMAJOR(ik))
      {
         pd->MB0 = ((mb0+ku-1)/ku)*ku;
         if (!ATL_AMM_KRUNTIME(ik))
            pd->amm_b0 = ATL_AMM_KERN_K1[ik];
      }
      else if (!ATL_AMM_KRUNTIME(ik) || mb0 != (mb0/ku)*ku ||
               mb0 < ATL_AMM_KBMINs[ik])
      {
         pd->amm_b0 = ATL_AMM_KERN_K1[ik];
         pd->MB0 = mb0;
      }
      else
         pd->MB0 = mb0;
      pd->nmu0 = (mb0+mu-1)/mu;
   }
   pd->mb0 = mb0;
   pd->nmblks = nmblks;
   pd->nxblks = nnblks * nmblks;
   if (P > pd->nxblks)
      P = pd->nxblks;
   pd->mu = mu;
   pd->nu = nu;
   pd->ku = ATL_AMM_KUs[ik];
   #ifdef TCPLX
      if (TA == AtlasConjTrans)
         pd->a2blk = ATL_AMM_AC2BLK_an[ik];
      else if (TA == AtlasConj)
         pd->a2blk = ATL_AMM_AH2BLK_an[ik];
      else
   #endif
   pd->a2blk = (TA == AtlasNoTrans) ?
               ATL_AMM_AT2BLK_an[ik] : ATL_AMM_A2BLK_an[ik];
   pd->b2blk = ATL_AMM_B2BLK_a1[ik];
/*
 * beta != 0, because then trsm simply zeros X and returns
 */
   if (SCALAR_IS_ONE(beta))
      pd->blk2c = ATL_AMM_BLK2C_a1_b1[ik];
   else if (SCALAR_IS_NONE(beta))
      pd->blk2c = ATL_AMM_BLK2C_a1_bn[ik];
   else
      pd->blk2c = ATL_AMM_BLK2C_a1_bX[ik];
// printf("ik=%d, mb=%d(%d), nb=%d(%d), MB0=%d\n", ik, pd->mb, pd->mb0, pd->nb, pd->nbf, pd->MB0);
   return(P);
}

ablk2cmat_t Mjoin(PATL,tGetSyammInfo)
(
   amminfo_t *out,
   const int P,          /* scale you want to use */
   enum ATLAS_TRANS TA,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alpha,
   const SCALAR beta
)
{
   ablk2cmat_t dblk2cmat;
   int ik = ATL_MAXIDX;
   int nb, k;
   if (K < ATL_AMM_MAXKB)
   {
      for (ik=0; ik <= ATL_MAXIDX && ATL_AMM_KBs[ik] < K; ik++);
   }
   nb = ATL_AMM_MBs[ik];
   #ifdef ATL_AMM_66IDX
      while (ik > ATL_AMM_66IDX)
      {
         k = N / nb;
         k = ((k-1)*k)>>1;
         if (k >= P)
            break;
         nb = ATL_AMM_MBs[--ik];
      }
   #endif
   out->IDX = ik;
   nb = Mmax(nb, ATL_AMM_NBs[ik]);
   if (nb > N)
      nb = N;
   out->mu = ATL_AMM_MUs[ik];
   out->nu = ATL_AMM_NUs[ik];
   out->ku = ATL_AMM_KUs[ik];
   k = ATL_lcm(out->mu, out->nu);
   nb = (nb > k) ? (nb/k)*k : k;
   out->nb = out->mb = nb;
   out->kb = ATL_AMM_KBs[ik];
/*   printf("tGetSyAMM, nb=%d, kb=%d, ik=%d\n", nb, out->kb, ik); */
   out->kbmin = ATL_AMM_KBMINs[ik];
   out->flag = ATL_AMM_KFLAG[ik];
   out->amm_b0 = ATL_AMM_KERN_b0[ik];
   out->amm_b1 = ATL_AMM_KERN_b1[ik];
   out->amm_bn = ATL_AMM_KERN_bn[ik];
   out->amm_k1_b0 = ATL_AMM_KERN_K1[ik];
   out->amm_k1_b1 = ATL_AMM_KERN_K1_b1[ik];
   out->amm_k1_bn = ATL_AMM_KERN_K1_bn[ik];
   if (TA == AtlasNoTrans)
   {
      out->a2blk = ATL_AMM_AT2BLK_a1[ik];
      out->b2blk =  ATL_AMM_BT2BLK_a1[ik];
   }
   else
   {
      out->a2blk = ATL_AMM_A2BLK_a1[ik];
      out->b2blk =  ATL_AMM_B2BLK_a1[ik];
   }
   if (SCALAR_IS_ONE(alpha))
   {
      dblk2cmat = ATL_AMM_BLK2C_a1_b0[ik];
      if (SCALAR_IS_ONE(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_a1_b1[ik];
      else if (SCALAR_IS_NONE(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_a1_bn[ik];
      else if (SCALAR_IS_ZERO(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_a1_b0[ik];
      else
         out->Cblk2cm = ATL_AMM_BLK2C_a1_bX[ik];
   }
   else if (SCALAR_IS_NONE(alpha))
   {
      dblk2cmat = ATL_AMM_BLK2C_an_b0[ik];
      if (SCALAR_IS_ONE(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_an_b1[ik];
      else if (SCALAR_IS_NONE(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_an_bn[ik];
      else if (SCALAR_IS_ZERO(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_an_b0[ik];
      else
         out->Cblk2cm = ATL_AMM_BLK2C_an_bX[ik];
   }
   else  /* alpha = X */
   {
      dblk2cmat = ATL_AMM_BLK2C_aX_b0[ik];
      if (SCALAR_IS_ONE(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_aX_b1[ik];
      else if (SCALAR_IS_NONE(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_aX_bn[ik];
      else if (SCALAR_IS_ZERO(beta))
         out->Cblk2cm = ATL_AMM_BLK2C_aX_b0[ik];
      else
         out->Cblk2cm = ATL_AMM_BLK2C_aX_bX[ik];
   }
   return(dblk2cmat);
}
/*
 * returns cblk2c_b0, cblk2c_b1 is in structure
 */
ablk2cmat_t Mjoin(PATL,tGetSyammInfo_K)
(
   amminfo_t *out,
   const int P,          /* scale you want to use */
   enum ATLAS_TRANS TA,
   ATL_CSZT N,
   ATL_CSZT K
)
{
   ablk2cmat_t dblk2cmat;
   int ik = ATL_AMM_NCASES-1;
   int mb, nb, k, mu, nu;

   if (K < ATL_AMM_MAXKB)
      for (ik=0; ik < ATL_AMM_NCASES-1 && ATL_AMM_KBs[ik] < K; ik++);
   out->IDX = ik;
   mu = out->mu = ATL_AMM_MUs[ik];
   nu = out->nu = ATL_AMM_NUs[ik];
   out->ku = ATL_AMM_KUs[ik];
   out->mb = ((N+mu-1)/mu)*mu;
   out->nb = ((N+nu-1)/nu)*nu;
   out->kb = ATL_AMM_KBs[ik];
/*  printf("tGetSyAMM_K, mb=%d, nb=%d, kb=%d, ik=%d\n", mb, nb, out->kb, ik); */
   out->kbmin = ATL_AMM_KBMINs[ik];
   out->flag = ATL_AMM_KFLAG[ik];
   out->amm_b0 = ATL_AMM_KERN_b0[ik];
   out->amm_b1 = ATL_AMM_KERN_b1[ik];
   out->amm_bn = ATL_AMM_KERN_bn[ik];
   out->amm_k1_b0 = ATL_AMM_KERN_K1[ik];
   out->amm_k1_b1 = ATL_AMM_KERN_K1_b1[ik];
   out->amm_k1_bn = ATL_AMM_KERN_K1_bn[ik];
   if (TA == AtlasNoTrans)
   {
      out->a2blk = ATL_AMM_AT2BLK_a1[ik];
      out->b2blk =  ATL_AMM_BT2BLK_a1[ik];
   }
   else
   {
      out->a2blk = ATL_AMM_A2BLK_a1[ik];
      out->b2blk =  ATL_AMM_B2BLK_a1[ik];
   }
   out->Cblk2cm = ATL_AMM_BLK2C_a1_b1[ik];
   return(ATL_AMM_BLK2C_a1_b0[ik]);
}
#endif

