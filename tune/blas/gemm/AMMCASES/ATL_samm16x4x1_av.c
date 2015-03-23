#include "atlas_misc.h"
#define ATL_NoFakePF
#include "atlas_prefetch.h"

#ifndef ATL_CSZT
   #define ATL_CSZT const size_t
#endif
void ATL_USERMM
(
   ATL_CSZT nmus,
   ATL_CSZT nnus,
   ATL_CSZT K,
   const TYPE *pA,    /* 4*KB*nmus-length access-major array of A */
   const TYPE *pB,    /* 1*KB*nnus-length access-major array of B */
   TYPE *pC,          /* 4*1*nnus*nmus-length access-major array of C */
   const TYPE *pAn,   /* next block of A */
   const TYPE *pBn,   /* next block of B */
   const TYPE *pCn    /* next block of C */
)
/*
 * Very basic 4x4 KU=1 AltiVec kernel
 */
{
   const TYPE *pB0 = pB, *pA0 = pA;
   const TYPE *pfA, *pfB;
   size_t incPF, i, j, k;
   vector float vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3;
   vector float  vC00, vC10, vC20, vC30, vC01, vC11, vC21, vC31,
                 vC02, vC12, vC22, vC32, vC03, vC13, vC23, vC33;
   #ifndef ATL_NoIEEE /* turn on java/ieee mode */
         const vector int izero   =  VECTOR_INITI(0,0,0,0);
         vec_mtvscr(izero);
   #endif

  if (pAn != pA)
   {
      pfA = pAn;
      incPF = (nmus*4*K) / (nmus * nnus);
      pfB = (pBn != pB) ? pBn : pCn;

   }
   else if (pCn != pC)
   {
      pfA = pCn;
      incPF = (nmus*4*nnus*1) / (nmus * nnus);
      pfB = pBn;
   }
   else if (pBn != pB)
   {
      pfA = pBn;
      pfB = pBn +((nmus*4*nnus*1)>>1);
      incPF = (K*nnus*1*sizeof(TYPE)) / ((nmus * nnus)<<1);
   }
   else
   {
      pfA = pA + nmus*4*(K>>1);
      incPF = (nmus*4*K) / (nmus * nnus);
   }
   for (i=0; i < nmus; i++)
   {
      for (j=0; j < nnus; j++)
      {
/*
 *       Peel K=0 iteration to prefetch
 */
         vC33 = vec_xor(vC33, vC33);
         vB3 = vec_ld(0, pB);
         ATL_pfl1W(pC);
         ATL_pfl1W(pC+32);
         vB0 = vec_splat(vB3, 0);
         vB1 = vec_splat(vB3, 1);
         vB2 = vec_splat(vB3, 2);
         vB3 = vec_splat(vB3, 3);

         vA0 = vec_ld(0, pA);
         vA1 = vec_ld(0, pA+4);
         vA2 = vec_ld(0, pA+8);
         vA3 = vec_ld(0, pA+12);

         vC00 = vec_madd(vA0, vB0, vC33);
         ATL_pfl1R(pfA);
         vC10 = vec_madd(vA1, vB0, vC33);
         vC20 = vec_madd(vA2, vB0, vC33);
         vC30 = vec_madd(vA3, vB0, vC33);
         vC01 = vec_madd(vA0, vB1, vC33);
         vC11 = vec_madd(vA1, vB1, vC33);
         vC21 = vec_madd(vA2, vB1, vC33);
         vC31 = vec_madd(vA3, vB1, vC33);
            vB1 = vec_ld(0, pB+4);
         vC02 = vec_madd(vA0, vB2, vC33);
            vB0 = vec_splat(vB1, 0);
         vC12 = vec_madd(vA1, vB2, vC33);
         vC22 = vec_madd(vA2, vB2, vC33);
         vC32 = vec_madd(vA3, vB2, vC33);
            vB2 = vec_splat(vB1, 2);
         vC03 = vec_madd(vA0, vB3, vC33);
            vA0 = vec_ld(0, pA+16);
         vC13 = vec_madd(vA1, vB3, vC33);
            vA1 = vec_ld(0, pA+20);
         vC23 = vec_madd(vA2, vB3, vC33);
            vA2 = vec_ld(0, pA+24);
         vC33 = vec_madd(vA3, vB3, vC33);
            vB3 = vec_splat(vB1, 3);
         pA += 32;
         pB += 8;
/*
 *       Handle remaining K its with rolled loop (compiler can unroll easily)
 */
      #if KB != 0
         for (k=1; k < KB; k++)
      #else
         for (k=1; k < K; k++)
      #endif
         {
               vA3 = vec_ld(0, pA-4);
            vC00 = vec_madd(vA0, vB0, vC00);
            vC10 = vec_madd(vA1, vB0, vC10);
            vC20 = vec_madd(vA2, vB0, vC20);
            vC30 = vec_madd(vA3, vB0, vC30);
               vB1 = vec_splat(vB1, 1);
            vC01 = vec_madd(vA0, vB1, vC01);
            vC11 = vec_madd(vA1, vB1, vC11);
            vC21 = vec_madd(vA2, vB1, vC21);
            vC31 = vec_madd(vA3, vB1, vC31);
               vB1 = vec_ld(0, pB);
            vC02 = vec_madd(vA0, vB2, vC02);
               vB0 = vec_splat(vB1, 0);
            vC12 = vec_madd(vA1, vB2, vC12);
            vC22 = vec_madd(vA2, vB2, vC22);
            vC32 = vec_madd(vA3, vB2, vC32);
               vB2 = vec_splat(vB1, 2);
            vC03 = vec_madd(vA0, vB3, vC03);
            vC13 = vec_madd(vA1, vB3, vC13);
               vA0 = vec_ld(0, pA);
               vA1 = vec_ld(0, pA+4);
            vC23 = vec_madd(vA2, vB3, vC23);
            vC33 = vec_madd(vA3, vB3, vC33);
               vA2 = vec_ld(0, pA+8);
               vB3 = vec_splat(vB1, 3);
            pA += 16;
            pB += 4;
         }
         #ifdef BETA0
           vec_st(vC00, 0, pC);
           vec_st(vC10, 0, pC+4);
           vec_st(vC20, 0, pC+8);
           vec_st(vC30, 0, pC+12);
           vec_st(vC01, 0, pC+16);
           vec_st(vC11, 0, pC+20);
           vec_st(vC21, 0, pC+24);
           vec_st(vC31, 0, pC+28);
           vec_st(vC02, 0, pC+32);
           vec_st(vC12, 0, pC+36);
           vec_st(vC22, 0, pC+40);
           vec_st(vC32, 0, pC+44);
           vec_st(vC03, 0, pC+48);
           vec_st(vC13, 0, pC+52);
           vec_st(vC23, 0, pC+56);
           vec_st(vC33, 0, pC+60);
         #else
            #ifdef BETAN1
               #define VEC_ADD vec_sub
            #else
               #define VEC_ADD vec_add
            #endif
            vA0 = vec_ld(0, pC);
            vA1 = vec_ld(0, pC+4);
            vA2 = vec_ld(0, pC+8);
            vA3 = vec_ld(0, pC+12);
            vC00 = VEC_ADD(vC00, vA0);
            vC10 = VEC_ADD(vC10, vA1);
            vC20 = VEC_ADD(vC20, vA2);
            vC30 = VEC_ADD(vC30, vA3);
            vec_st(vC00, 0, pC);
            vec_st(vC10, 0, pC+4);
            vec_st(vC20, 0, pC+8);
            vec_st(vC30, 0, pC+12);
            vB0 = vec_ld(0, pC+16);
            vB1 = vec_ld(0, pC+20);
            vB2 = vec_ld(0, pC+24);
            vB3 = vec_ld(0, pC+28);
            vC01 = VEC_ADD(vC01, vB0);
            vC11 = VEC_ADD(vC11, vB1);
            vC21 = VEC_ADD(vC21, vB2);
            vC31 = VEC_ADD(vC31, vB3);
            vec_st(vC01, 0, pC+16);
            vec_st(vC11, 0, pC+20);
            vec_st(vC21, 0, pC+24);
            vec_st(vC31, 0, pC+28);
            vA0 = vec_ld(0, pC+32);
            vA1 = vec_ld(0, pC+36);
            vA2 = vec_ld(0, pC+40);
            vA3 = vec_ld(0, pC+44);
            vC02 = VEC_ADD(vC02, vA0);
            vC12 = VEC_ADD(vC12, vA1);
            vC22 = VEC_ADD(vC22, vA2);
            vC32 = VEC_ADD(vC32, vA3);
            vec_st(vC02, 0, pC+32);
            vec_st(vC12, 0, pC+36);
            vec_st(vC22, 0, pC+40);
            vec_st(vC32, 0, pC+44);
            vB0 = vec_ld(0, pC+48);
            vB1 = vec_ld(0, pC+52);
            vB2 = vec_ld(0, pC+56);
            vB3 = vec_ld(0, pC+60);
            vC03 = VEC_ADD(vC03, vB0);
            vC13 = VEC_ADD(vC13, vB1);
            vC23 = VEC_ADD(vC23, vB2);
            vC33 = VEC_ADD(vC33, vB3);
            vec_st(vC03, 0, pC+48);
            vec_st(vC13, 0, pC+52);
            vec_st(vC23, 0, pC+56);
            vec_st(vC33, 0, pC+60);
         #endif
         pA = pA0;
         pC += 16*4;  /* MU * NU */
         pB -= 4;
      }
      pA0 += 16*K;     /* MU*K */
      pA = pA0;
      pB = pB0;
   }
}
