#include <stdio.h>
#include <stdlib.h>
int main(int nargs, char **args)
{
   double *x, *y, *z, ans[8];
   void *vp;
   void do_vmacc(double *z, double *x, double *y);
   int i;

   vp = malloc(8*3*sizeof(double) + 64);
   x = (double*) ( 64 + ((((size_t)(vp))>>6)<<6) );
   y = x + 8;
   z = y + 8;
   x[0] = 1.0; x[1] = 2.0; x[2] = 4.0; x[3] = 8.0;
   x[4] = 16.0; x[5] = 32.0; x[6] = 64.0; x[7] = 128.0;
   y[0] = -4.0; y[1] = -8.0; y[2] = -16.0; y[3] = -32.0;
   y[4] = -64.0; y[5] = -128.0; y[6] = -256.0; y[7] = -512.0;
   z[0] = 0.0; z[1] = 1.0; z[2] = -1.0; z[3] = -4.0;
   z[4] = 4.0; z[5] = -8.0; z[6] = 8.0; z[7] = -16.0;
   for (i=0; i < 8; i++)
      ans[i] = z[i] + x[i] * y[i];
   do_vmacc(z, x, y);   /* z += x * y */
   if (z[0] != ans[0] || z[1] != ans[1] || z[2] != ans[2] || z[3] != ans[3] ||
       z[4] != ans[4] || z[5] != ans[5] || z[6] != ans[6] || z[7] != ans[7])
   {
      fprintf(stderr, "wanted={%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f}\n",
              ans[0],ans[1],ans[2],ans[3],ans[4],ans[5],ans[6],ans[7]);
      fprintf(stderr, "got   ={%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f}\n",
              z[0],z[1],z[2],z[3],z[4],z[5],z[6],z[7]);
      printf("FAILURE\n");
      free(vp);
      exit(1);
   }
   printf("SUCCESS\n");
   free(vp);
   return(0);
}
