#include "atlas_tlapack.h"
#include "atlas_pca.h"

void Mjoin(PATL,DoWorkLASWP)(ATL_LAUNCHSTRUCT_t *lp, void *vp)
{
   ATL_thread_t *tp=vp;
   ATL_TLASWP_N_t *swpp=((ATL_TLASWP_N_t*)lp->opstruct)+tp->rank;
   const int *piv=swpp->ipiv, *ipiv;
   ATL_CINT nr=swpp->nr, K1=swpp->K1, K2=swpp->K2, inci=swpp->inci;
   ATL_INT nblks=swpp->nblks;
   TYPE *A = swpp->A;
   size_t lda = (swpp->lda)SHIFT;
   const int n = K2 - K1;
   size_t incA = lda << 5;
   int i, ip, i1, i2, KeepOn;
   register int h;
   TYPE *a0, *a1;
   #ifdef TCPLX
      register TYPE r0, r1;
   #else
      register TYPE r;
   #endif

   if (K2 < K1) return;
   if (inci < 0)
   {
      piv -= (K2-1) * inci;
      i1 = K2-1;
      i2 = K1;
   }
   else
   {
      piv += K1*inci;
      i1 = K1;
      i2 = K2-1;
   }

   if (nblks)
   {
      do
      {
         ipiv = piv;
         i = i1;
         do
         {
            ip = *ipiv; ipiv += inci;
            if (ip != i)
            {
               a0 = A + (i SHIFT);
               a1 = A + (ip SHIFT);
               for (h=32; h; h--)
               {
                  #ifdef TCPLX
                     r0 = *a0;
                     r1 = a0[1];
                     *a0 = *a1;
                     a0[1] = a1[1];
                     *a1 = r0;
                     a1[1] = r1;
                  #else
                     r = *a0;
                     *a0 = *a1;
                     *a1 = r;
                  #endif
                  a0 += lda;
                  a1 += lda;
               }
            }
            if (inci > 0) KeepOn = (++i <= i2);
            else KeepOn = (--i >= i2);
         }
         while(KeepOn);
         A += incA;
      }
      while(--nblks);
   }
   if (nr)
   {
      ipiv = piv;
      i = i1;
      do
      {
         ip = *ipiv; ipiv += inci;
         if (ip != i)
         {
            a0 = A + (i SHIFT);
            a1 = A + (ip SHIFT);
            for (h=nr; h; h--)
            {
               #ifdef TCPLX
                  r0 = *a0;
                  r1 = a0[1];
                  *a0 = *a1;
                  a0[1] = a1[1];
                  *a1 = r0;
                  a1[1] = r1;
               #else
                  r = *a0;
                  *a0 = *a1;
                  *a1 = r;
               #endif
               a0 += lda;
               a1 += lda;
            }
         }
         if (inci > 0) KeepOn = (++i <= i2);
         else KeepOn = (--i >= i2);
      }
      while(KeepOn);
   }
}

void Mjoin(PATL,tlaswp)(const int N, TYPE *A, const int lda, const int K1,
                        const int K2, const int *ipiv, const int inci)
{
   ATL_TLASWP_N_t swps[ATL_NTHREADS];
   ATL_INT i, p = ATL_NTHREADS, nblks, nlblks, neblks, nr, n;
   TYPE *ac;
   if (N < 128)
   {
      Mjoin(PATL,laswp)(N, A, lda, K1, K2, ipiv, inci);
      return;
   }
   nblks = N >> 5;
   nr = N - (nblks<<5);
   if (nblks < p)
   {
      p = nblks;
      nlblks = 1;
      neblks = 0;
   }
   else
   {
      nlblks = nblks / p;
      neblks = nblks - nlblks*p;
   }
   ac = A;
   for (i=0; i < p; i++)
   {
      size_t n;
      n = swps[i].nblks = (i < neblks) ? nlblks+1 : nlblks;
      swps[i].nr = (i == neblks) ? nr : 0;
      n = (n<<5)+swps[i].nr;
      swps[i].A = ac;
      swps[i].K1 = K1;
      swps[i].K2 = K2;
      swps[i].ipiv = ipiv;
      swps[i].inci = inci;
      swps[i].lda = lda;
      ac += (((size_t)lda)SHIFT)*n;
   }
   ATL_goparallel(p, Mjoin(PATL,DoWorkLASWP), swps, NULL);
}
