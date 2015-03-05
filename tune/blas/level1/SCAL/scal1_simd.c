#include "atlas_misc.h"
#include "atlas_simd.h"
#include "atlas_prefetch.h"
#define XDIST 128
void ATL_USCAL(const int N, const SCALAR alpha, TYPE *X, const int incX)
{
   ATL_VTYPE valp;
   ATL_CSZT vlmskb = (ATL_VLENb-1);
   int i;
   int nr, nv;

   ATL_vbcast(valp, &alpha);
/*
 * Peel loop until aligned to VLEN
 */
   for (i=0; i < N && ((size_t)(X)&vlmskb); i++)
      *X++ *= alpha;

   nr = N-i;
   nv = nr & (~(ATL_VLEN-1));
   if (nv)
   {
      nr -= nv;
      nv += i;
      for (; i < nv; i += ATL_VLEN, X += ATL_VLEN)
      {
         ATL_VTYPE x0;
         ATL_vld(x0, X);
         ATL_vmul(x0, x0, valp);
            ATL_pfl1R(X+XDIST);
         ATL_vst(X, x0);
      }
   }
   for (; i < N; i++)
      *X++ *= alpha;
}
