#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),amm_rankK.h))
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_rkkblk.h))
#include Mstr(Mjoin(Mjoin(atlas_,UPR),amm_rkkflag.h))

int Mjoin(PATL,GetRankKInfo)
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
   const int ik = K-3;
   int appAl;  /* 0:A, 1:B */
   ATL_assert(K > 2 && K <= ATL_MAXK_RKK);
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
   out->amm_b0 = ATL_AMM_KERN_RKK[ik];
   out->amm_b1 = ATL_AMM_KERN_RKK_b1[ik];
   out->amm_bn = ATL_AMM_KERN_RKK_bn[ik];
/*
 * Apply alpha to smallest matrix, and use alpha/beta to pick copy routines
 */
   if (SCALAR_IS_ONE(alpha))
   {
      appAl = 0;
      #ifdef TCPLX
         if (TA == AtlasNoTrans)
            out->a2blk = ATL_RKK_AT2BLK_a1[ik];
         else if (TA == AtlasTrans)
            out->a2blk = ATL_RKK_A2BLK_a1[ik];
         else if (TA == AtlasConjTrans)
            out->a2blk = ATL_RKK_AC2BLK_a1[ik];
         else
            out->a2blk = ATL_RKK_AH2BLK_a1[ik];
         if (TB == AtlasNoTrans)
             out->b2blk = ATL_RKK_B2BLK_a1[ik];
         else if (TB == AtlasTrans)
             out->b2blk = ATL_RKK_BT2BLK_a1[ik];
         else if (TB == AtlasConjTrans)
             out->b2blk = ATL_RKK_BH2BLK_a1[ik];
         else
             out->b2blk = ATL_RKK_BC2BLK_a1[ik];
      #else
         out->a2blk = (TA == AtlasNoTrans) ?
            ATL_RKK_AT2BLK_a1[ik]:ATL_RKK_A2BLK_a1[ik];
         out->b2blk = (TB == AtlasNoTrans) ?
            ATL_RKK_B2BLK_a1[ik]:ATL_RKK_BT2BLK_a1[ik];
      #endif
      if (SCALAR_IS_ONE(beta))
         out->Cblk2cm = ATL_RKK_BLK2C_a1_b1[ik];
      else if (SCALAR_IS_ZERO(beta))
         out->Cblk2cm = ATL_RKK_BLK2C_a1_b0[ik];
      else if (SCALAR_IS_NONE(beta))
         out->Cblk2cm = ATL_RKK_BLK2C_a1_bn[ik];
      else
         out->Cblk2cm = ATL_RKK_BLK2C_a1_bX[ik];
   }
   else  /* alpha is not one */
   {
      appAl = (M >= N) ? 1:0;
      if (SCALAR_IS_ONE(beta))
         out->Cblk2cm = ATL_RKK_BLK2C_a1_b1[ik];
      else if (SCALAR_IS_ZERO(beta))
         out->Cblk2cm = ATL_RKK_BLK2C_a1_b0[ik];
      else if (SCALAR_IS_NONE(beta))
         out->Cblk2cm = ATL_RKK_BLK2C_a1_bn[ik];
      else
         out->Cblk2cm = ATL_RKK_BLK2C_a1_bX[ik];
      if (!appAl)  /* apply to alpha to A */
      {
         #ifdef TCPLX
            if (TB == AtlasNoTrans)
                out->b2blk = ATL_RKK_B2BLK_a1[ik];
            else if (TB == AtlasTrans)
                out->b2blk = ATL_RKK_BT2BLK_a1[ik];
            else if (TB == AtlasConjTrans)
                out->b2blk = ATL_RKK_BH2BLK_a1[ik];
            else
                out->b2blk = ATL_RKK_BC2BLK_a1[ik];
            if (SCALAR_IS_NONE(alpha))
            {
               if (TA == AtlasNoTrans)
                  out->a2blk = ATL_RKK_AT2BLK_an[ik];
               else if (TA == AtlasTrans)
                  out->a2blk = ATL_RKK_A2BLK_an[ik];
               else if (TA == AtlasConjTrans)
                  out->a2blk = ATL_RKK_AC2BLK_an[ik];
               else
                  out->a2blk = ATL_RKK_AH2BLK_an[ik];
            }
            else
            {
               if (TA == AtlasNoTrans)
                  out->a2blk = ATL_RKK_AT2BLK_aX[ik];
               else if (TA == AtlasTrans)
                  out->a2blk = ATL_RKK_A2BLK_aX[ik];
               else if (TA == AtlasConjTrans)
                  out->a2blk = ATL_RKK_AC2BLK_aX[ik];
               else
                  out->a2blk = ATL_RKK_AH2BLK_aX[ik];
            }
         #else
            if (SCALAR_IS_NONE(alpha))
               out->a2blk = (TA == AtlasNoTrans) ?
                            ATL_RKK_AT2BLK_an[ik] : ATL_RKK_A2BLK_an[ik];
            else
               out->a2blk = (TA == AtlasNoTrans) ?
                            ATL_RKK_AT2BLK_aX[ik] : ATL_RKK_A2BLK_aX[ik];
            out->b2blk = (TB == AtlasNoTrans) ?
                         ATL_RKK_B2BLK_a1[ik] : ATL_RKK_BT2BLK_a1[ik];
         #endif
      }
      else /* apply alpha to B */
      {
         #ifdef TCPLX
            if (TA == AtlasNoTrans)
               out->a2blk = ATL_RKK_AT2BLK_a1[ik];
            else if (TA == AtlasTrans)
               out->a2blk = ATL_RKK_A2BLK_a1[ik];
            else if (TA == AtlasConjTrans)
               out->a2blk = ATL_RKK_AC2BLK_a1[ik];
            else
               out->a2blk = ATL_RKK_AH2BLK_a1[ik];
            if (SCALAR_IS_NONE(alpha))
            {
               if (TB == AtlasNoTrans)
                   out->b2blk = ATL_RKK_B2BLK_an[ik];
               else if (TB == AtlasTrans)
                   out->b2blk = ATL_RKK_BT2BLK_an[ik];
               else if (TB == AtlasConjTrans)
                   out->b2blk = ATL_RKK_BH2BLK_an[ik];
               else
                   out->b2blk = ATL_RKK_BC2BLK_an[ik];
            }
            else
            {
               if (TB == AtlasNoTrans)
                   out->b2blk = ATL_RKK_B2BLK_aX[ik];
               else if (TB == AtlasTrans)
                   out->b2blk = ATL_RKK_BT2BLK_aX[ik];
               else if (TB == AtlasConjTrans)
                   out->b2blk = ATL_RKK_BH2BLK_aX[ik];
               else
                   out->b2blk = ATL_RKK_BC2BLK_aX[ik];
            }
         #else
            out->a2blk = (TA == AtlasNoTrans) ?
                         ATL_RKK_AT2BLK_a1[ik] : ATL_RKK_A2BLK_a1[ik];
            if (SCALAR_IS_NONE(alpha))
               out->b2blk = (TB == AtlasNoTrans) ?
                            ATL_RKK_B2BLK_an[ik] : ATL_RKK_BT2BLK_an[ik];
            else
               out->b2blk = (TB == AtlasNoTrans) ?
                            ATL_RKK_B2BLK_aX[ik] : ATL_RKK_BT2BLK_aX[ik];
         #endif
      }
   }
   return(appAl);
}

