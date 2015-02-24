#include "atlas_asm.h"
#
# x86 assembler for:
# void do_vmacc(z, x, y)
# RETURNS: z[0:7] += x[0:7] *  y[0:7];
#
.text
.globl ATL_asmdecor(do_vmacc)
ATL_asmdecor(do_vmacc):
   #ifdef ATL_GAS_x8664
      vmovapd 0(%rdi), %zmm0                  /* load z */
      vmovapd 0(%rsi), %zmm1                  /* load x */
      vmovapd 0(%rdx), %zmm2                  /* load y */
      vfmadd231pd %zmm2, %zmm1, %zmm0         /* z += x*y */
      vmovapd %zmm0, (%rdi)                   /* store z */
   #else
      movl    12(%esp), %eax
      vmovapd (%eax), %zmm2                   /* load y */
      movl    8(%esp), %eax
      vmovapd  (%eax), %zmm1                  /* load x */
      movl    4(%esp), %eax
      vmovapd (%eax), %zmm0                   /* load z */
      vfmadd231pd %zmm2, %zmm1, %zmm0         /* z += x*y */
      vmovapd %zmm0, (%eax)                   /* store z */
   #endif
   ret
