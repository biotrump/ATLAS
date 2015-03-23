#include "atlas_misc.h"
#include "atlas_simd.h"
#if defined(__GNUC__) || \
    (defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999))
   #define ATL_SINLINE static inline
#else
   #define ATL_SINLINE static
#endif
#if ATL_VLEN > 1
ATL_SINLINE void ATL_rk4n4_uvec
(
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   const SCALAR alpha,
   const TYPE *A0,
   ATL_CSZT lda,
   const TYPE *B0,
   ATL_CSZT ldb,
   TYPE *C0,
   ATL_CSZT ldc
)
{
   ATL_VTYPE b00, b10, b20, b30, b01, b11, b21, b31;
   ATL_VTYPE b02, b12, b22, b32, b03, b13, b23, b33;
   const TYPE *A1=A0+lda, *A2=A1+lda, *A3=A2+lda;
   TYPE *C1=C0+ldc, *C2=C1+ldc, *C3=C2+ldc;
   register int i;
   if (TB == AtlasNoTrans)
   {
      const TYPE *B1=B0+ldb, *B2=B1+ldb, *B3=B2+ldb;
      ATL_vbcast(b00, B0);
      ATL_vbcast(b01, B1);
      ATL_vbcast(b02, B2);
      ATL_vbcast(b03, B3);
      ATL_vbcast(b10, B0+1);
      ATL_vbcast(b11, B1+1);
      ATL_vbcast(b12, B2+1);
      ATL_vbcast(b13, B3+1);
      ATL_vbcast(b20, B0+2);
      ATL_vbcast(b21, B1+2);
      ATL_vbcast(b22, B2+2);
      ATL_vbcast(b23, B3+2);
      ATL_vbcast(b30, B0+3);
      ATL_vbcast(b31, B1+3);
      ATL_vbcast(b32, B2+3);
      ATL_vbcast(b33, B3+3);
   }
   else /* TB == AtlasTrans, B NxK */
   {
      const TYPE *B1=B0+ldb, *B2=B1+ldb, *B3=B2+ldb;
      ATL_vbcast(b00, B0);
      ATL_vbcast(b10, B1);
      ATL_vbcast(b20, B2);
      ATL_vbcast(b30, B3);
      ATL_vbcast(b01, B0+1);
      ATL_vbcast(b11, B1+1);
      ATL_vbcast(b21, B2+1);
      ATL_vbcast(b31, B3+1);
      ATL_vbcast(b02, B0+2);
      ATL_vbcast(b12, B1+2);
      ATL_vbcast(b22, B2+2);
      ATL_vbcast(b32, B3+2);
      ATL_vbcast(b03, B0+3);
      ATL_vbcast(b13, B1+3);
      ATL_vbcast(b23, B2+3);
      ATL_vbcast(b33, B3+3);
   }
   if (alpha != 1.0)
   {
      ATL_VTYPE al;
      ATL_vbcast(al, &alpha);
      ATL_vmul(b00, b00, al);
      ATL_vmul(b10, b10, al);
      ATL_vmul(b20, b20, al);
      ATL_vmul(b30, b30, al);
      ATL_vmul(b01, b01, al);
      ATL_vmul(b11, b11, al);
      ATL_vmul(b21, b21, al);
      ATL_vmul(b31, b31, al);
      ATL_vmul(b02, b02, al);
      ATL_vmul(b12, b12, al);
      ATL_vmul(b22, b22, al);
      ATL_vmul(b32, b32, al);
      ATL_vmul(b03, b03, al);
      ATL_vmul(b13, b13, al);
      ATL_vmul(b23, b23, al);
      ATL_vmul(b33, b33, al);
   }
   for (i=0; i < M; i += ATL_VLEN)
   {
      ATL_VTYPE c0, c1, c2, c3, a0;
      ATL_vuld(c0, C0+i);
      ATL_vuld(c1, C1+i);
      ATL_vuld(c2, C2+i);
      ATL_vuld(c3, C3+i);

      ATL_vuld(a0, A0+i);
      ATL_vmac(c0, b00, a0);
      ATL_vmac(c1, b01, a0);
      ATL_vmac(c2, b02, a0);
      ATL_vmac(c3, b03, a0);
      ATL_vuld(a0, A1+i);
      ATL_vmac(c0, b10, a0);
      ATL_vmac(c1, b11, a0);
      ATL_vmac(c2, b12, a0);
      ATL_vmac(c3, b13, a0);
      ATL_vuld(a0, A2+i);
      ATL_vmac(c0, b20, a0);
      ATL_vmac(c1, b21, a0);
      ATL_vmac(c2, b22, a0);
      ATL_vmac(c3, b23, a0);
      ATL_vuld(a0, A3+i);
      ATL_vmac(c0, b30, a0);
      ATL_vmac(c1, b31, a0);
      ATL_vmac(c2, b32, a0);
      ATL_vmac(c3, b33, a0);

      ATL_vust(C0+i, c0);
      ATL_vust(C1+i, c1);
      ATL_vust(C2+i, c2);
      ATL_vust(C3+i, c3);
   }
}
ATL_SINLINE void ATL_rk4n4_vec
(
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   const SCALAR alpha,
   const TYPE *A0,
   ATL_CSZT lda,
   const TYPE *B0,
   ATL_CSZT ldb,
   TYPE *C0,
   ATL_CSZT ldc
)
{
   ATL_VTYPE b00, b10, b20, b30, b01, b11, b21, b31;
   ATL_VTYPE b02, b12, b22, b32, b03, b13, b23, b33;
   const TYPE *A1=A0+lda, *A2=A1+lda, *A3=A2+lda;
   TYPE *C1=C0+ldc, *C2=C1+ldc, *C3=C2+ldc;
   register int i;
   if (TB == AtlasNoTrans)
   {
      const TYPE *B1=B0+ldb, *B2=B1+ldb, *B3=B2+ldb;
      ATL_vbcast(b00, B0);
      ATL_vbcast(b01, B1);
      ATL_vbcast(b02, B2);
      ATL_vbcast(b03, B3);
      ATL_vbcast(b10, B0+1);
      ATL_vbcast(b11, B1+1);
      ATL_vbcast(b12, B2+1);
      ATL_vbcast(b13, B3+1);
      ATL_vbcast(b20, B0+2);
      ATL_vbcast(b21, B1+2);
      ATL_vbcast(b22, B2+2);
      ATL_vbcast(b23, B3+2);
      ATL_vbcast(b30, B0+3);
      ATL_vbcast(b31, B1+3);
      ATL_vbcast(b32, B2+3);
      ATL_vbcast(b33, B3+3);
   }
   else /* TB == AtlasTrans, B NxK */
   {
      const TYPE *B1=B0+ldb, *B2=B1+ldb, *B3=B2+ldb;
      ATL_vbcast(b00, B0);
      ATL_vbcast(b10, B1);
      ATL_vbcast(b20, B2);
      ATL_vbcast(b30, B3);
      ATL_vbcast(b01, B0+1);
      ATL_vbcast(b11, B1+1);
      ATL_vbcast(b21, B2+1);
      ATL_vbcast(b31, B3+1);
      ATL_vbcast(b02, B0+2);
      ATL_vbcast(b12, B1+2);
      ATL_vbcast(b22, B2+2);
      ATL_vbcast(b32, B3+2);
      ATL_vbcast(b03, B0+3);
      ATL_vbcast(b13, B1+3);
      ATL_vbcast(b23, B2+3);
      ATL_vbcast(b33, B3+3);
   }
   if (alpha != 1.0)
   {
      ATL_VTYPE al;
      ATL_vbcast(al, &alpha);
      ATL_vmul(b00, b00, al);
      ATL_vmul(b10, b10, al);
      ATL_vmul(b20, b20, al);
      ATL_vmul(b30, b30, al);
      ATL_vmul(b01, b01, al);
      ATL_vmul(b11, b11, al);
      ATL_vmul(b21, b21, al);
      ATL_vmul(b31, b31, al);
      ATL_vmul(b02, b02, al);
      ATL_vmul(b12, b12, al);
      ATL_vmul(b22, b22, al);
      ATL_vmul(b32, b32, al);
      ATL_vmul(b03, b03, al);
      ATL_vmul(b13, b13, al);
      ATL_vmul(b23, b23, al);
      ATL_vmul(b33, b33, al);
   }
   for (i=0; i < M; i += ATL_VLEN)
   {
      ATL_VTYPE c0, c1, c2, c3, a0;
      ATL_vld(c0, C0+i);
      ATL_vld(c1, C1+i);
      ATL_vld(c2, C2+i);
      ATL_vld(c3, C3+i);

      ATL_vld(a0, A0+i);
      ATL_vmac(c0, b00, a0);
      ATL_vmac(c1, b01, a0);
      ATL_vmac(c2, b02, a0);
      ATL_vmac(c3, b03, a0);
      ATL_vld(a0, A1+i);
      ATL_vmac(c0, b10, a0);
      ATL_vmac(c1, b11, a0);
      ATL_vmac(c2, b12, a0);
      ATL_vmac(c3, b13, a0);
      ATL_vld(a0, A2+i);
      ATL_vmac(c0, b20, a0);
      ATL_vmac(c1, b21, a0);
      ATL_vmac(c2, b22, a0);
      ATL_vmac(c3, b23, a0);
      ATL_vld(a0, A3+i);
      ATL_vmac(c0, b30, a0);
      ATL_vmac(c1, b31, a0);
      ATL_vmac(c2, b32, a0);
      ATL_vmac(c3, b33, a0);

      ATL_vst(C0+i, c0);
      ATL_vst(C1+i, c1);
      ATL_vst(C2+i, c2);
      ATL_vst(C3+i, c3);
   }
}
#endif

ATL_SINLINE void ATL_rk4n4
(
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   const SCALAR alpha,
   const TYPE *A0,
   ATL_CSZT lda,
   const TYPE *B0,
   ATL_CSZT ldb,
   TYPE *C0,
   ATL_CSZT ldc
)

/*
 * This special-case code used for LU, only called when:
 *   TA=AtlasNoTrans, N == K == 4, beta == 1.0
 */
{
   TYPE b00, b10, b20, b30, b01, b11, b21, b31;
   TYPE b02, b12, b22, b32, b03, b13, b23, b33;
   const TYPE *A1=A0+lda, *A2=A1+lda, *A3=A2+lda;
   TYPE *C1=C0+ldc, *C2=C1+ldc, *C3=C2+ldc;
   register int i;

   if (TB == AtlasNoTrans)
   {
      const TYPE *B1=B0+ldb, *B2=B1+ldb, *B3=B2+ldb;
      b00 = alpha * *B0;
      b01 = alpha * *B1;
      b02 = alpha * *B2;
      b03 = alpha * *B3;
      b10 = alpha * B0[1];
      b11 = alpha * B1[1];
      b12 = alpha * B2[1];
      b13 = alpha * B3[1];
      b20 = alpha * B0[2];
      b21 = alpha * B1[2];
      b22 = alpha * B2[2];
      b23 = alpha * B3[2];
      b30 = alpha * B0[3];
      b31 = alpha * B1[3];
      b32 = alpha * B2[3];
      b33 = alpha * B3[3];
   }
   else  /* B == AtlasTrans; B is NxK */
   {
      const TYPE *B1=B0+ldb, *B2=B1+ldb, *B3=B2+ldb;
      b00 = alpha * *B0;
      b10 = alpha * *B1;
      b20 = alpha * *B2;
      b30 = alpha * *B3;
      b01 = alpha * B0[1];
      b11 = alpha * B1[1];
      b21 = alpha * B2[1];
      b31 = alpha * B3[1];
      b02 = alpha * B0[2];
      b12 = alpha * B1[2];
      b22 = alpha * B2[2];
      b32 = alpha * B3[2];
      b03 = alpha * B0[3];
      b13 = alpha * B1[3];
      b23 = alpha * B2[3];
      b33 = alpha * B3[3];
   }
   for (i=0; i < M; i++)
   {
      const register TYPE a0=A0[i], a1=A1[i], a2=A2[i], a3=A3[i];
      TYPE register c0=C0[i], c1=C1[i], c2=C2[i], c3=C3[i];
      c0 += b00*a0;
      c1 += b01*a0;
      c2 += b02*a0;
      c3 += b03*a0;
      c0 += b10*a1;
      c1 += b11*a1;
      c2 += b12*a1;
      c3 += b13*a1;
      c0 += b20*a2;
      c1 += b21*a2;
      c2 += b22*a2;
      c3 += b23*a2;
      c0 += b30*a3;
      c1 += b31*a3;
      c2 += b32*a3;
      c3 += b33*a3;
      C0[i] = c0;
      C1[i] = c1;
      C2[i] = c2;
      C3[i] = c3;
   }
}

int Mjoin(PATL,rk4n4)
(
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   const SCALAR alpha,
   const TYPE *A0,
   ATL_CSZT lda,
   const TYPE *B0,
   ATL_CSZT ldb,
   TYPE *C0,
   ATL_CSZT ldc
)
{
   #if ATL_VLEN < 2
      ATL_rk4n4(TB, M, alpha, A0, lda, B0, ldb, C0, ldc);
   #else
      const int vmod = (ATL_VLEN-1);
      int Mp=0, Mv, Mr;  /* peel, vector, remainder */
      int ALLALIGN=0;
      if (!((ldc&vmod) | (lda&vmod)))
      {
         const size_t VMOD = ((size_t)(ATL_VLEN-1));
         size_t a0, a1;
         int gap0, gap;
         a0 = (size_t) C0;
         a1 = (size_t) A0;
         a0 = ATL_DivBySize(a0);
         a1 = ATL_DivBySize(a1);
         gap0 = a0 & VMOD;
         gap  = a1 & VMOD;
         if (gap0 == gap)
         {
            ALLALIGN=1;
            if (gap)
            {
               Mp = ATL_VLEN - gap;
               Mp = Mmin(Mp, M);
            }
         }
      }
      Mr = M - Mp;
      Mv = Mr & ~vmod;
      Mr -= Mv;

      if (ALLALIGN)
      {
         if (Mp)
         {
            ATL_rk4n4(TB, Mp, alpha, A0, lda, B0, ldb, C0, ldc);
            A0 += Mp; C0 += Mp;
         }
         if (Mv)
         {
            ATL_rk4n4_vec(TB, Mv, alpha, A0, lda, B0, ldb, C0, ldc);
            A0 += Mv; C0 += Mv;
         }
      }
      else if (Mv)
      {
         ATL_rk4n4_uvec(TB, Mv, alpha, A0, lda, B0, ldb, C0, ldc);
         A0 += Mv; C0 += Mv;
      }
      if (Mr)
         ATL_rk4n4(TB, Mr, alpha, A0, lda, B0, ldb, C0, ldc);
   #endif
   return(0);
}
