#include "atlas_misc.h"
#include "atlas_simd.h"
#include "atlas_prefetch.h"
TYPE ATL_UDOT(const int N, const TYPE *X, const int incX,
             const TYPE *Y, const int incY)
{
   register TYPE dot=ATL_rzero;
   ATL_CSZT vlmskb = (ATL_VLENb-1);
   TYPE sum[ATL_VLEN];
   int i, nr, nv;

   for (i=0; i < N && ((size_t)(Y)&vlmskb); i++)
      dot += *X++ * *Y++;
   nr = N-i;
   nv = nr >> (ATL_VLSH+2);
   if (nv)
   {
      ATL_VTYPE vd0, vd1, vd2, vd3;
      int k;

      ATL_vzero(vd0);
      ATL_vzero(vd1);
      ATL_vzero(vd2);
      ATL_vzero(vd3);
      nv <<= (ATL_VLSH+2);
      nr -= nv;
      nv += i;
      if (0 &&((size_t)(X)&vlmskb) == 0)  /* Can we use aligned lds for X & Y? */
      {
         for (; i < nv; i += 4*ATL_VLEN, X += 4*ATL_VLEN, Y += 4*ATL_VLEN)
         {
            ATL_VTYPE x0, y0;
            ATL_vld(x0, X);
            ATL_vld(y0, Y);
            ATL_vmac(vd0, x0, y0);
            ATL_vld(x0, X+ATL_VLEN);
            ATL_vld(y0, Y+ATL_VLEN);
            ATL_vmac(vd1, x0, y0);
            ATL_vld(x0, X+2*ATL_VLEN);
            ATL_vld(y0, Y+2*ATL_VLEN);
            ATL_vmac(vd2, x0, y0);
            ATL_vld(x0, X+3*ATL_VLEN);
            ATL_vld(y0, Y+3*ATL_VLEN);
            ATL_vmac(vd3, x0, y0);
               ATL_pfl1R(Y+192);
         }
      }
      else
      {
         for (; i < nv; i += 4*ATL_VLEN, X += 4*ATL_VLEN, Y += 4*ATL_VLEN)
         {
            ATL_VTYPE x0, y0;
            ATL_vuld(x0, X);
            ATL_vld(y0, Y);
            ATL_vmac(vd0, x0, y0);
            ATL_vuld(x0, X+ATL_VLEN);
            ATL_vld(y0, Y+ATL_VLEN);
            ATL_vmac(vd1, x0, y0);
            ATL_vuld(x0, X+2*ATL_VLEN);
            ATL_vld(y0, Y+2*ATL_VLEN);
            ATL_vmac(vd2, x0, y0);
            ATL_vuld(x0, X+3*ATL_VLEN);
            ATL_vld(y0, Y+3*ATL_VLEN);
            ATL_vmac(vd3, x0, y0);
               ATL_pfl1R(Y+192);
//               ATL_pfl1R(X+64);
         }
      }
      ATL_vadd(vd0, vd0, vd1);
      ATL_vadd(vd2, vd2, vd3);
      ATL_vadd(vd0, vd0, vd2);
      ATL_vust(sum, vd0);
      for (k=0; k < ATL_VLEN; k++)
         dot += sum[k];
   }
   for (; i < N; i++)
      dot += *X++ * *Y++;
   return(dot);
}
