#include <stdlib.h>
#include "atlas_misc.h"

void ATL_ipivinit(ATL_UINT N, int *ipiv)
{
   int i;
   for (i=0; i < N; i++)
   {
      const unsigned int n = N-i;
      ipiv[i] = i + rand()%n;
   }
}

int *ATL_getIpiv(ATL_UINT N, ATL_UINT seed)
{
   int *ip;
   srand(seed);
   ip = malloc(N*sizeof(int));
   ATL_assert(ip);
   ATL_ipivinit(N, ip);
   return(ip);
}
