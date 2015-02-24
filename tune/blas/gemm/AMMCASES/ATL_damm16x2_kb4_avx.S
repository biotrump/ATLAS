#include "atlas_asm.h"
#define m0      %ymm0
#define rA0     %ymm1
#define rA1     %ymm2
#define rA2     %ymm3
#define rA3     %ymm4
#define rB0     %ymm5
#define rB1     %ymm6
#define rb1     %ymm7
#define rC00    %ymm8
#define rC10    %ymm9
#define rC20    %ymm10
#define rC30    %ymm11
#define rC01    %ymm12
#define rC11    %ymm13
#define rC21    %ymm14
#define rC31    %ymm15
#ifdef BETA1
   #define ADDC(r, p) vaddpd p, r, r
#elif defined(BETA0)
   #define ADDC(r, p)
#else
   #define ADDC(r, p) vsubpd p, r, r
#endif

#define NMU %rdi
#define NNU %rsi
#define pA  %rcx
#define pB  %rdx
#define pC  %r9
#define pB0 %r8
#define NNU0 %rax
/*
                    rdi      rsi    rdx        rcx         r8        r9
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                  8(%rsp)    16(%rsp)     24(%rsp)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */
.text
.global ATL_asmdecor(ATL_USERMM)
ALIGN16
ATL_asmdecor(ATL_USERMM):
   mov pB0, pB
   mov NNU, NNU0
   MLOOP:
      vbroadcastsd (pB), rB0
      vmovapd (pA), rA0
      NLOOP:
/*
 *       K == 1
*/
         vmulpd rA0, rB0, rC00
         ADDC(rC00, (pC))
         vmovapd 32(pA), rA1
         vmulpd rA1, rB0, rC10
         ADDC(rC10, 32(pC))
         vmovapd 64(pA), rA2
         vmulpd rA2, rB0, rC20
         ADDC(rC20, 64(pC))
         vmovapd 96(pA), rA3
         vmulpd rA3, rB0, rC30
         ADDC(rC30, 96(pC))
         vbroadcastsd 8(pB), rB1
         vmulpd rA0, rB1, rC01
         ADDC(rC01, 128(pC))
            vbroadcastsd 16(pB), rB0
         vmulpd rA1, rB1, rC11
         ADDC(rC11, 160(pC))
            vmovapd 128(pA), rA0
         vmulpd rA2, rB1, rC21
         ADDC(rC21, 192(pC))
            vmovapd 160(pA), rA1
         vmulpd rA3, rB1, rC31
         ADDC(rC31, 224(pC))
            vmovapd 192(pA), rA2
/*
 *       K == 2
 */
         vmulpd rA0, rB0, m0
         vaddpd m0, rC00, rC00
            vmovapd 224(pA), rA3
         vmulpd rA1, rB0, m0
         vaddpd m0, rC10, rC10
            vbroadcastsd 24(pB), rB1
         vmulpd rA2, rB0, m0
         vaddpd m0, rC20, rC20
         vmulpd rA3, rB0, m0
         vaddpd m0, rC30, rC30
            vbroadcastsd 32(pB), rB0
         vmulpd rA0, rB1, m0
         vaddpd m0, rC01, rC01
            vmovapd 256(pA), rA0
         vmulpd rA1, rB1, m0
         vaddpd m0, rC11, rC11
            vmovapd 288(pA), rA1
         vmulpd rA2, rB1, m0
         vaddpd m0, rC21, rC21
            vmovapd 320(pA), rA2
         vmulpd rA3, rB1, m0
         vaddpd m0, rC31, rC31
            vmovapd 352(pA), rA3
/*
 *       K == 3
 */
         vmulpd rA0, rB0, m0
         vaddpd m0, rC00, rC00
            vbroadcastsd 40(pB), rB1
         vmulpd rA1, rB0, m0
         vaddpd m0, rC10, rC10
            vbroadcastsd 56(pB), rb1
         vmulpd rA2, rB0, m0
         vaddpd m0, rC20, rC20
         vmulpd rA3, rB0, m0
         vaddpd m0, rC30, rC30
            vbroadcastsd 48(pB), rB0
         vmulpd rA0, rB1, m0
         vaddpd m0, rC01, rC01
            vmovapd 384(pA), rA0
         vmulpd rA1, rB1, m0
         vaddpd m0, rC11, rC11
            vmovapd 416(pA), rA1
         vmulpd rA2, rB1, m0
         vaddpd m0, rC21, rC21
            vmovapd 448(pA), rA2
         vmulpd rA3, rB1, m0
         vaddpd m0, rC31, rC31
            vmovapd 480(pA), rA3
/*
 *       K == 4
 */
         vmulpd rA0, rB0, m0
         vaddpd m0, rC00, rC00
         vmovapd rC00, (pC)
         vmulpd rA1, rB0, m0
         vaddpd m0, rC10, rC10
         vmovapd rC10, 32(pC)
         vmulpd rA2, rB0, m0
         vaddpd m0, rC20, rC20
         vmovapd rC20, 64(pC)
         vmulpd rA3, rB0, m0
         vaddpd m0, rC30, rC30
         vmovapd rC30, 96(pC)
            vbroadcastsd 64(pB), rB0
         add $64, pB
         vmulpd rA0, rb1, m0
         vaddpd m0, rC01, rC01
         vmovapd rC01, 128(pC)
            vmovapd (pA), rA0
         vmulpd rA1, rb1, m0
         vaddpd m0, rC11, rC11
         vmovapd rC11, 160(pC)
         vmulpd rA2, rb1, m0
         vaddpd m0, rC21, rC21
         vmovapd rC21, 192(pC)
         vmulpd rA3, rb1, m0
         vaddpd m0, rC31, rC31
         vmovapd rC31, 224(pC)
         add $256, pC
      sub $1, NNU
      jnz NLOOP

      add $512, pA
      mov pB0, pB
      mov NNU0, NNU
   sub $1, NMU
   jnz MLOOP
/* DONE: */
   ret
