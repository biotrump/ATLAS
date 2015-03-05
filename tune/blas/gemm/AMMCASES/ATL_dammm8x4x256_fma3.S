#include "atlas_asm.h"

#define rB0     %ymm0
#define rB1     %ymm1
#define rB2     %ymm2
#define rB3     %ymm3
#define rA0     %ymm4
#define rA1     %ymm5
#define ra0     %ymm6
#define ra1     %ymm7
#define rC00    %ymm8
#define rC10    %ymm9
#define rC01    %ymm10
#define rC11    %ymm11
#define rC02    %ymm12
#define rC12    %ymm13
#define rC03    %ymm14
#define rC13    %ymm15
/*
 * Prioritize original registers for inner-loop operations, but inc regs
 * can be anything w/o changing opcode size, so use new regs for those
 */
#define KK      %rdx  /* API reg */
#define pA      %rcx  /* API reg */
#define pB      %rax  /* comes in as r9 */
#define r24     %r9   /* set after mov r9 to pC () */
/*
 * Then N-loop variables much less important, so use any orig regs left
 */
#define pA0     %r8   /* set after mov r8 to pB (rax) */
#define pC      %rsi  /* set after mov rsi to nnu () */
#define nnu     %r10  /* comes in as rsi */
#define pfA     %rbx
#define pfB     %rbp
#define incPF   %r12
#define KK0     %rdi
/*
 * We could give a rat's ass about what registers used in outer (M-) loop
 */
#define nmu     %r11  /* comes in as rdi */
#define incAm   %r13
#define nnu0    %r14
#define pB0     %r15
/*
                    rdi      rsi    rdx        rcx         r8        r9
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                  8(%rsp)    16(%rsp)     24(%rsp)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */
#define PFBDIST 128
#define PFADIST 128
#define prefA(m_) prefetcht2 m_
#define prefB(m_) prefetcht2 m_
#define prefC(m_) prefetcht0 m_
#define FMAC vfmadd231pd   /* FMAC m256/r256, rs1, rd */
#if defined(BETAN) || defined(BETAn)
   #define BETAN1
#endif
#ifdef BETAN1
   #define VCOP vsubpd
#else
   #define VCOP vaddpd
#endif
#define movapd movaps
.text
.global ATL_asmdecor(ATL_USERMM)
ALIGN16
ATL_asmdecor(ATL_USERMM):
/*
 * Save callee-saved iregs
 */
   movq %rbp, -8(%rsp)
   movq %rbx, -16(%rsp)
   movq %r12, -24(%rsp)
   movq %r13, -32(%rsp)
   movq %r14, -40(%rsp)
   movq %r15, -48(%rsp)
/*
 * Load paramaters
 */
   mov %rdi, nmu
   mov %rsi, nnu
   mov %r8, pB
   mov %r9, pC
   mov nnu, nnu0
   movq 8(%rsp), pfB      /* pfB = pAn */
   movq 16(%rsp), pfA     /* pf = pBn */
   cmp pfB, pB
   CMOVE pfA, pfB
   mov KK, KK0
   sub $-128, pfB
   sub $-128, pfA
   sub $-128, pA
   sub $-128, pB
   sub $-128, pC
   mov $24, r24
   mov pA, pA0
   mov pB, pB0
   mov $8*4*8, incPF
/*
 * incAm = 8*sizeof*K = 8*8*K = 64*K
 */
   mov KK, incAm
   shl $6, incAm           /* incAm = 32*K */

   ALIGN16
   MLOOP:
      NLOOP:
/*
 *       First iteration peeled to handle init of rC
 */
         vxorps rC00, rC00, rC00
         vmovaps rC00, rC10
         vxorps rC01, rC01, rC01
         vmovaps rC00, rC11
         vxorps rC02, rC02, rC02
         vmovaps rC00, rC12
         vxorps rC03, rC03, rC03
         vmovaps rC00, rC13

         vmovaps -128(pA), rA0

         #if 0
         vmovaps -128(pB), rb0          /* rb0 = {b3, b2, b1, b0} */
         vunpcklpd rb0, rb0, rb1        /* rb1 = {b2, b2, b0, b0} */
         vperf2f128 $0, rb1, rb1, rB0   /* rB0 = {b0, b0, b0, b0} */
         vunpckhpd rb0, rb0, rb0        /* rb0 = {b3, b3  b1, b1} */
         vperf2f128 $0, rb0, rb0, rB1   /* rB1 = {b1, b1, b1, b1} */

         vperf2f128 0x11, rb1, rb1, rB0 /* rB0 = {b2, b2, b2, b2} */
         vperf2f128 0x11, rb0, rb0, rB1 /* rB1 = {b3, b3, b3, b3} */
         #else
         vbroadcastsd -128(pB), rB0
         vbroadcastsd -120(pB), rB1
         vbroadcastsd -112(pB), rB2
         #endif

         vmovaps -96(pA), rA1
         #if KB > 1
            vmovaps -64(pA), ra0
            vmovaps -32(pA), ra1
         #endif
/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         ALIGN16
         #if KB > 0
            vbroadcastsd -104(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 1
               vbroadcastsd -96(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 1
               vbroadcastsd -88(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 1
               vbroadcastsd -80(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 2
               vmovaps 0(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 2
               vmovaps 32(pA), rA1
            #endif
         #endif
         #if KB > 1
            vbroadcastsd -72(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 2
               vbroadcastsd -64(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 2
               vbroadcastsd -56(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 2
               vbroadcastsd -48(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 3
               vmovaps 64(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 3
               vmovaps 96(pA), ra1
            #endif
         #endif
         #if KB > 2
            vbroadcastsd -40(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 3
               vbroadcastsd -32(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 3
               vbroadcastsd -24(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 3
               vbroadcastsd -16(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 4
               vmovaps 128(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 4
               vmovaps 160(pA), rA1
            #endif
         #endif
         #if KB > 3
            vbroadcastsd -8(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 4
               vbroadcastsd 0(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 4
               vbroadcastsd 8(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 4
               vbroadcastsd 16(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 5
               vmovaps 192(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 5
               vmovaps 224(pA), ra1
            #endif
         #endif
         #if KB > 4
            vbroadcastsd 24(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 5
               vbroadcastsd 32(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 5
               vbroadcastsd 40(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 5
               vbroadcastsd 48(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 6
               vmovaps 256(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 6
               vmovaps 288(pA), rA1
            #endif
         #endif
         #if KB > 5
            vbroadcastsd 56(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 6
               vbroadcastsd 64(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 6
               vbroadcastsd 72(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 6
               vbroadcastsd 80(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 7
               vmovaps 320(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 7
               vmovaps 352(pA), ra1
            #endif
         #endif
         #if KB > 6
            vbroadcastsd 88(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 7
               vbroadcastsd 96(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 7
               vbroadcastsd 104(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 7
               vbroadcastsd 112(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 8
               vmovaps 384(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 8
               vmovaps 416(pA), rA1
            #endif
         #endif
         #if KB > 7
            vbroadcastsd 120(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 8
               vbroadcastsd 128(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 8
               vbroadcastsd 136(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 8
               vbroadcastsd 144(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 9
               vmovaps 448(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 9
               vmovaps 480(pA), ra1
            #endif
         #endif
         #if KB > 8
            vbroadcastsd 152(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 9
               vbroadcastsd 160(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 9
               vbroadcastsd 168(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 9
               vbroadcastsd 176(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 10
               vmovaps 512(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 10
               vmovaps 544(pA), rA1
            #endif
         #endif
         #if KB > 9
            vbroadcastsd 184(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 10
               vbroadcastsd 192(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 10
               vbroadcastsd 200(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 10
               vbroadcastsd 208(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 11
               vmovaps 576(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 11
               vmovaps 608(pA), ra1
            #endif
         #endif
         #if KB > 10
            vbroadcastsd 216(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 11
               vbroadcastsd 224(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 11
               vbroadcastsd 232(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 11
               vbroadcastsd 240(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 12
               vmovaps 640(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 12
               vmovaps 672(pA), rA1
            #endif
         #endif
         #if KB > 11
            vbroadcastsd 248(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 12
               vbroadcastsd 256(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 12
               vbroadcastsd 264(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 12
               vbroadcastsd 272(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 13
               vmovaps 704(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 13
               vmovaps 736(pA), ra1
            #endif
         #endif
         #if KB > 12
            vbroadcastsd 280(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 13
               vbroadcastsd 288(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 13
               vbroadcastsd 296(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 13
               vbroadcastsd 304(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 14
               vmovaps 768(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 14
               vmovaps 800(pA), rA1
            #endif
         #endif
         #if KB > 13
            vbroadcastsd 312(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 14
               vbroadcastsd 320(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 14
               vbroadcastsd 328(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 14
               vbroadcastsd 336(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 15
               vmovaps 832(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 15
               vmovaps 864(pA), ra1
            #endif
         #endif
         #if KB > 14
            vbroadcastsd 344(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 15
               vbroadcastsd 352(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 15
               vbroadcastsd 360(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 15
               vbroadcastsd 368(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 16
               vmovaps 896(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 16
               vmovaps 928(pA), rA1
            #endif
         #endif
         #if KB > 15
            vbroadcastsd 376(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 16
               vbroadcastsd 384(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 16
               vbroadcastsd 392(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 16
               vbroadcastsd 400(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 17
               vmovaps 960(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 17
               vmovaps 992(pA), ra1
            #endif
         #endif
         #if KB > 16
            vbroadcastsd 408(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 17
               vbroadcastsd 416(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 17
               vbroadcastsd 424(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 17
               vbroadcastsd 432(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 18
               vmovaps 1024(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 18
               vmovaps 1056(pA), rA1
            #endif
         #endif
         #if KB > 17
            vbroadcastsd 440(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 18
               vbroadcastsd 448(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 18
               vbroadcastsd 456(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 18
               vbroadcastsd 464(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 19
               vmovaps 1088(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 19
               vmovaps 1120(pA), ra1
            #endif
         #endif
         #if KB > 18
            vbroadcastsd 472(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 19
               vbroadcastsd 480(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 19
               vbroadcastsd 488(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 19
               vbroadcastsd 496(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 20
               vmovaps 1152(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 20
               vmovaps 1184(pA), rA1
            #endif
         #endif
         #if KB > 19
            vbroadcastsd 504(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 20
               vbroadcastsd 512(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 20
               vbroadcastsd 520(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 20
               vbroadcastsd 528(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 21
               vmovaps 1216(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 21
               vmovaps 1248(pA), ra1
            #endif
         #endif
         #if KB > 20
            vbroadcastsd 536(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 21
               vbroadcastsd 544(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 21
               vbroadcastsd 552(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 21
               vbroadcastsd 560(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 22
               vmovaps 1280(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 22
               vmovaps 1312(pA), rA1
            #endif
         #endif
         #if KB > 21
            vbroadcastsd 568(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 22
               vbroadcastsd 576(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 22
               vbroadcastsd 584(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 22
               vbroadcastsd 592(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 23
               vmovaps 1344(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 23
               vmovaps 1376(pA), ra1
            #endif
         #endif
         #if KB > 22
            vbroadcastsd 600(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 23
               vbroadcastsd 608(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 23
               vbroadcastsd 616(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 23
               vbroadcastsd 624(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 24
               vmovaps 1408(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 24
               vmovaps 1440(pA), rA1
            #endif
         #endif
         #if KB > 23
            vbroadcastsd 632(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 24
               vbroadcastsd 640(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 24
               vbroadcastsd 648(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 24
               vbroadcastsd 656(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 25
               vmovaps 1472(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 25
               vmovaps 1504(pA), ra1
            #endif
         #endif
         #if KB > 24
            vbroadcastsd 664(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 25
               vbroadcastsd 672(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 25
               vbroadcastsd 680(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 25
               vbroadcastsd 688(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 26
               vmovaps 1536(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 26
               vmovaps 1568(pA), rA1
            #endif
         #endif
         #if KB > 25
            vbroadcastsd 696(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 26
               vbroadcastsd 704(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 26
               vbroadcastsd 712(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 26
               vbroadcastsd 720(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 27
               vmovaps 1600(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 27
               vmovaps 1632(pA), ra1
            #endif
         #endif
         #if KB > 26
            vbroadcastsd 728(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 27
               vbroadcastsd 736(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 27
               vbroadcastsd 744(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 27
               vbroadcastsd 752(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 28
               vmovaps 1664(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 28
               vmovaps 1696(pA), rA1
            #endif
         #endif
         #if KB > 27
            vbroadcastsd 760(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 28
               vbroadcastsd 768(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 28
               vbroadcastsd 776(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 28
               vbroadcastsd 784(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 29
               vmovaps 1728(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 29
               vmovaps 1760(pA), ra1
            #endif
         #endif
         #if KB > 28
            vbroadcastsd 792(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 29
               vbroadcastsd 800(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 29
               vbroadcastsd 808(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 29
               vbroadcastsd 816(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 30
               vmovaps 1792(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 30
               vmovaps 1824(pA), rA1
            #endif
         #endif
         #if KB > 29
            vbroadcastsd 824(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 30
               vbroadcastsd 832(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 30
               vbroadcastsd 840(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 30
               vbroadcastsd 848(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 31
               vmovaps 1856(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 31
               vmovaps 1888(pA), ra1
            #endif
         #endif
         #if KB > 30
            vbroadcastsd 856(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 31
               vbroadcastsd 864(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 31
               vbroadcastsd 872(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 31
               vbroadcastsd 880(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 32
               vmovaps 1920(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 32
               vmovaps 1952(pA), rA1
            #endif
         #endif
         #if KB > 31
            vbroadcastsd 888(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 32
               vbroadcastsd 896(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 32
               vbroadcastsd 904(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 32
               vbroadcastsd 912(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 33
               vmovaps 1984(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 33
               vmovaps 2016(pA), ra1
            #endif
         #endif
         #if KB > 32
            vbroadcastsd 920(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 33
               vbroadcastsd 928(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 33
               vbroadcastsd 936(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 33
               vbroadcastsd 944(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 34
               vmovaps 2048(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 34
               vmovaps 2080(pA), rA1
            #endif
         #endif
         #if KB > 33
            vbroadcastsd 952(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 34
               vbroadcastsd 960(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 34
               vbroadcastsd 968(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 34
               vbroadcastsd 976(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 35
               vmovaps 2112(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 35
               vmovaps 2144(pA), ra1
            #endif
         #endif
         #if KB > 34
            vbroadcastsd 984(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 35
               vbroadcastsd 992(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 35
               vbroadcastsd 1000(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 35
               vbroadcastsd 1008(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 36
               vmovaps 2176(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 36
               vmovaps 2208(pA), rA1
            #endif
         #endif
         #if KB > 35
            vbroadcastsd 1016(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 36
               vbroadcastsd 1024(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 36
               vbroadcastsd 1032(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 36
               vbroadcastsd 1040(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 37
               vmovaps 2240(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 37
               vmovaps 2272(pA), ra1
            #endif
         #endif
         #if KB > 36
            vbroadcastsd 1048(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 37
               vbroadcastsd 1056(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 37
               vbroadcastsd 1064(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 37
               vbroadcastsd 1072(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 38
               vmovaps 2304(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 38
               vmovaps 2336(pA), rA1
            #endif
         #endif
         #if KB > 37
            vbroadcastsd 1080(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 38
               vbroadcastsd 1088(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 38
               vbroadcastsd 1096(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 38
               vbroadcastsd 1104(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 39
               vmovaps 2368(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 39
               vmovaps 2400(pA), ra1
            #endif
         #endif
         #if KB > 38
            vbroadcastsd 1112(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 39
               vbroadcastsd 1120(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 39
               vbroadcastsd 1128(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 39
               vbroadcastsd 1136(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 40
               vmovaps 2432(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 40
               vmovaps 2464(pA), rA1
            #endif
         #endif
         #if KB > 39
            vbroadcastsd 1144(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 40
               vbroadcastsd 1152(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 40
               vbroadcastsd 1160(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 40
               vbroadcastsd 1168(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 41
               vmovaps 2496(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 41
               vmovaps 2528(pA), ra1
            #endif
         #endif
         #if KB > 40
            vbroadcastsd 1176(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 41
               vbroadcastsd 1184(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 41
               vbroadcastsd 1192(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 41
               vbroadcastsd 1200(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 42
               vmovaps 2560(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 42
               vmovaps 2592(pA), rA1
            #endif
         #endif
         #if KB > 41
            vbroadcastsd 1208(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 42
               vbroadcastsd 1216(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 42
               vbroadcastsd 1224(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 42
               vbroadcastsd 1232(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 43
               vmovaps 2624(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 43
               vmovaps 2656(pA), ra1
            #endif
         #endif
         #if KB > 42
            vbroadcastsd 1240(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 43
               vbroadcastsd 1248(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 43
               vbroadcastsd 1256(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 43
               vbroadcastsd 1264(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 44
               vmovaps 2688(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 44
               vmovaps 2720(pA), rA1
            #endif
         #endif
         #if KB > 43
            vbroadcastsd 1272(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 44
               vbroadcastsd 1280(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 44
               vbroadcastsd 1288(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 44
               vbroadcastsd 1296(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 45
               vmovaps 2752(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 45
               vmovaps 2784(pA), ra1
            #endif
         #endif
         #if KB > 44
            vbroadcastsd 1304(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 45
               vbroadcastsd 1312(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 45
               vbroadcastsd 1320(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 45
               vbroadcastsd 1328(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 46
               vmovaps 2816(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 46
               vmovaps 2848(pA), rA1
            #endif
         #endif
         #if KB > 45
            vbroadcastsd 1336(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 46
               vbroadcastsd 1344(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 46
               vbroadcastsd 1352(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 46
               vbroadcastsd 1360(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 47
               vmovaps 2880(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 47
               vmovaps 2912(pA), ra1
            #endif
         #endif
         #if KB > 46
            vbroadcastsd 1368(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 47
               vbroadcastsd 1376(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 47
               vbroadcastsd 1384(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 47
               vbroadcastsd 1392(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 48
               vmovaps 2944(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 48
               vmovaps 2976(pA), rA1
            #endif
         #endif
         #if KB > 47
            vbroadcastsd 1400(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 48
               vbroadcastsd 1408(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 48
               vbroadcastsd 1416(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 48
               vbroadcastsd 1424(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 49
               vmovaps 3008(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 49
               vmovaps 3040(pA), ra1
            #endif
         #endif
         #if KB > 48
            vbroadcastsd 1432(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 49
               vbroadcastsd 1440(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 49
               vbroadcastsd 1448(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 49
               vbroadcastsd 1456(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 50
               vmovaps 3072(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 50
               vmovaps 3104(pA), rA1
            #endif
         #endif
         #if KB > 49
            vbroadcastsd 1464(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 50
               vbroadcastsd 1472(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 50
               vbroadcastsd 1480(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 50
               vbroadcastsd 1488(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 51
               vmovaps 3136(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 51
               vmovaps 3168(pA), ra1
            #endif
         #endif
         #if KB > 50
            vbroadcastsd 1496(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 51
               vbroadcastsd 1504(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 51
               vbroadcastsd 1512(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 51
               vbroadcastsd 1520(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 52
               vmovaps 3200(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 52
               vmovaps 3232(pA), rA1
            #endif
         #endif
         #if KB > 51
            vbroadcastsd 1528(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 52
               vbroadcastsd 1536(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 52
               vbroadcastsd 1544(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 52
               vbroadcastsd 1552(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 53
               vmovaps 3264(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 53
               vmovaps 3296(pA), ra1
            #endif
         #endif
         #if KB > 52
            vbroadcastsd 1560(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 53
               vbroadcastsd 1568(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 53
               vbroadcastsd 1576(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 53
               vbroadcastsd 1584(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 54
               vmovaps 3328(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 54
               vmovaps 3360(pA), rA1
            #endif
         #endif
         #if KB > 53
            vbroadcastsd 1592(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 54
               vbroadcastsd 1600(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 54
               vbroadcastsd 1608(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 54
               vbroadcastsd 1616(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 55
               vmovaps 3392(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 55
               vmovaps 3424(pA), ra1
            #endif
         #endif
         #if KB > 54
            vbroadcastsd 1624(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 55
               vbroadcastsd 1632(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 55
               vbroadcastsd 1640(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 55
               vbroadcastsd 1648(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 56
               vmovaps 3456(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 56
               vmovaps 3488(pA), rA1
            #endif
         #endif
         #if KB > 55
            vbroadcastsd 1656(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 56
               vbroadcastsd 1664(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 56
               vbroadcastsd 1672(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 56
               vbroadcastsd 1680(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 57
               vmovaps 3520(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 57
               vmovaps 3552(pA), ra1
            #endif
         #endif
         #if KB > 56
            vbroadcastsd 1688(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 57
               vbroadcastsd 1696(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 57
               vbroadcastsd 1704(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 57
               vbroadcastsd 1712(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 58
               vmovaps 3584(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 58
               vmovaps 3616(pA), rA1
            #endif
         #endif
         #if KB > 57
            vbroadcastsd 1720(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 58
               vbroadcastsd 1728(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 58
               vbroadcastsd 1736(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 58
               vbroadcastsd 1744(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 59
               vmovaps 3648(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 59
               vmovaps 3680(pA), ra1
            #endif
         #endif
         #if KB > 58
            vbroadcastsd 1752(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 59
               vbroadcastsd 1760(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 59
               vbroadcastsd 1768(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 59
               vbroadcastsd 1776(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 60
               vmovaps 3712(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 60
               vmovaps 3744(pA), rA1
            #endif
         #endif
         #if KB > 59
            vbroadcastsd 1784(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 60
               vbroadcastsd 1792(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 60
               vbroadcastsd 1800(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 60
               vbroadcastsd 1808(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 61
               vmovaps 3776(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 61
               vmovaps 3808(pA), ra1
            #endif
         #endif
         #if KB > 60
            vbroadcastsd 1816(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 61
               vbroadcastsd 1824(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 61
               vbroadcastsd 1832(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 61
               vbroadcastsd 1840(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 62
               vmovaps 3840(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 62
               vmovaps 3872(pA), rA1
            #endif
         #endif
         #if KB > 61
            vbroadcastsd 1848(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 62
               vbroadcastsd 1856(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 62
               vbroadcastsd 1864(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 62
               vbroadcastsd 1872(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 63
               vmovaps 3904(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 63
               vmovaps 3936(pA), ra1
            #endif
         #endif
         #if KB > 62
            vbroadcastsd 1880(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 63
               vbroadcastsd 1888(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 63
               vbroadcastsd 1896(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 63
               vbroadcastsd 1904(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 64
               vmovaps 3968(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 64
               vmovaps 4000(pA), rA1
            #endif
         #endif
         #if KB > 63
            vbroadcastsd 1912(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 64
               vbroadcastsd 1920(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 64
               vbroadcastsd 1928(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 64
               vbroadcastsd 1936(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 65
               vmovaps 4032(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 65
               vmovaps 4064(pA), ra1
            #endif
         #endif
         #if KB > 64
            vbroadcastsd 1944(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 65
               vbroadcastsd 1952(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 65
               vbroadcastsd 1960(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 65
               vbroadcastsd 1968(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 66
               vmovaps 4096(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 66
               vmovaps 4128(pA), rA1
            #endif
         #endif
         #if KB > 65
            vbroadcastsd 1976(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 66
               vbroadcastsd 1984(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 66
               vbroadcastsd 1992(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 66
               vbroadcastsd 2000(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 67
               vmovaps 4160(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 67
               vmovaps 4192(pA), ra1
            #endif
         #endif
         #if KB > 66
            vbroadcastsd 2008(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 67
               vbroadcastsd 2016(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 67
               vbroadcastsd 2024(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 67
               vbroadcastsd 2032(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 68
               vmovaps 4224(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 68
               vmovaps 4256(pA), rA1
            #endif
         #endif
         #if KB > 67
            vbroadcastsd 2040(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 68
               vbroadcastsd 2048(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 68
               vbroadcastsd 2056(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 68
               vbroadcastsd 2064(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 69
               vmovaps 4288(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 69
               vmovaps 4320(pA), ra1
            #endif
         #endif
         #if KB > 68
            vbroadcastsd 2072(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 69
               vbroadcastsd 2080(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 69
               vbroadcastsd 2088(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 69
               vbroadcastsd 2096(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 70
               vmovaps 4352(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 70
               vmovaps 4384(pA), rA1
            #endif
         #endif
         #if KB > 69
            vbroadcastsd 2104(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 70
               vbroadcastsd 2112(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 70
               vbroadcastsd 2120(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 70
               vbroadcastsd 2128(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 71
               vmovaps 4416(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 71
               vmovaps 4448(pA), ra1
            #endif
         #endif
         #if KB > 70
            vbroadcastsd 2136(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 71
               vbroadcastsd 2144(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 71
               vbroadcastsd 2152(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 71
               vbroadcastsd 2160(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 72
               vmovaps 4480(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 72
               vmovaps 4512(pA), rA1
            #endif
         #endif
         #if KB > 71
            vbroadcastsd 2168(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 72
               vbroadcastsd 2176(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 72
               vbroadcastsd 2184(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 72
               vbroadcastsd 2192(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 73
               vmovaps 4544(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 73
               vmovaps 4576(pA), ra1
            #endif
         #endif
         #if KB > 72
            vbroadcastsd 2200(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 73
               vbroadcastsd 2208(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 73
               vbroadcastsd 2216(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 73
               vbroadcastsd 2224(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 74
               vmovaps 4608(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 74
               vmovaps 4640(pA), rA1
            #endif
         #endif
         #if KB > 73
            vbroadcastsd 2232(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 74
               vbroadcastsd 2240(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 74
               vbroadcastsd 2248(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 74
               vbroadcastsd 2256(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 75
               vmovaps 4672(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 75
               vmovaps 4704(pA), ra1
            #endif
         #endif
         #if KB > 74
            vbroadcastsd 2264(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 75
               vbroadcastsd 2272(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 75
               vbroadcastsd 2280(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 75
               vbroadcastsd 2288(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 76
               vmovaps 4736(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 76
               vmovaps 4768(pA), rA1
            #endif
         #endif
         #if KB > 75
            vbroadcastsd 2296(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 76
               vbroadcastsd 2304(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 76
               vbroadcastsd 2312(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 76
               vbroadcastsd 2320(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 77
               vmovaps 4800(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 77
               vmovaps 4832(pA), ra1
            #endif
         #endif
         #if KB > 76
            vbroadcastsd 2328(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 77
               vbroadcastsd 2336(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 77
               vbroadcastsd 2344(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 77
               vbroadcastsd 2352(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 78
               vmovaps 4864(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 78
               vmovaps 4896(pA), rA1
            #endif
         #endif
         #if KB > 77
            vbroadcastsd 2360(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 78
               vbroadcastsd 2368(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 78
               vbroadcastsd 2376(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 78
               vbroadcastsd 2384(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 79
               vmovaps 4928(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 79
               vmovaps 4960(pA), ra1
            #endif
         #endif
         #if KB > 78
            vbroadcastsd 2392(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 79
               vbroadcastsd 2400(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 79
               vbroadcastsd 2408(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 79
               vbroadcastsd 2416(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 80
               vmovaps 4992(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 80
               vmovaps 5024(pA), rA1
            #endif
         #endif
         #if KB > 79
            vbroadcastsd 2424(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 80
               vbroadcastsd 2432(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 80
               vbroadcastsd 2440(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 80
               vbroadcastsd 2448(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 81
               vmovaps 5056(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 81
               vmovaps 5088(pA), ra1
            #endif
         #endif
         #if KB > 80
            vbroadcastsd 2456(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 81
               vbroadcastsd 2464(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 81
               vbroadcastsd 2472(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 81
               vbroadcastsd 2480(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 82
               vmovaps 5120(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 82
               vmovaps 5152(pA), rA1
            #endif
         #endif
         #if KB > 81
            vbroadcastsd 2488(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 82
               vbroadcastsd 2496(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 82
               vbroadcastsd 2504(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 82
               vbroadcastsd 2512(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 83
               vmovaps 5184(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 83
               vmovaps 5216(pA), ra1
            #endif
         #endif
         #if KB > 82
            vbroadcastsd 2520(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 83
               vbroadcastsd 2528(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 83
               vbroadcastsd 2536(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 83
               vbroadcastsd 2544(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 84
               vmovaps 5248(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 84
               vmovaps 5280(pA), rA1
            #endif
         #endif
         #if KB > 83
            vbroadcastsd 2552(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 84
               vbroadcastsd 2560(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 84
               vbroadcastsd 2568(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 84
               vbroadcastsd 2576(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 85
               vmovaps 5312(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 85
               vmovaps 5344(pA), ra1
            #endif
         #endif
         #if KB > 84
            vbroadcastsd 2584(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 85
               vbroadcastsd 2592(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 85
               vbroadcastsd 2600(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 85
               vbroadcastsd 2608(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 86
               vmovaps 5376(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 86
               vmovaps 5408(pA), rA1
            #endif
         #endif
         #if KB > 85
            vbroadcastsd 2616(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 86
               vbroadcastsd 2624(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 86
               vbroadcastsd 2632(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 86
               vbroadcastsd 2640(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 87
               vmovaps 5440(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 87
               vmovaps 5472(pA), ra1
            #endif
         #endif
         #if KB > 86
            vbroadcastsd 2648(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 87
               vbroadcastsd 2656(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 87
               vbroadcastsd 2664(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 87
               vbroadcastsd 2672(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 88
               vmovaps 5504(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 88
               vmovaps 5536(pA), rA1
            #endif
         #endif
         #if KB > 87
            vbroadcastsd 2680(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 88
               vbroadcastsd 2688(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 88
               vbroadcastsd 2696(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 88
               vbroadcastsd 2704(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 89
               vmovaps 5568(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 89
               vmovaps 5600(pA), ra1
            #endif
         #endif
         #if KB > 88
            vbroadcastsd 2712(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 89
               vbroadcastsd 2720(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 89
               vbroadcastsd 2728(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 89
               vbroadcastsd 2736(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 90
               vmovaps 5632(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 90
               vmovaps 5664(pA), rA1
            #endif
         #endif
         #if KB > 89
            vbroadcastsd 2744(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 90
               vbroadcastsd 2752(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 90
               vbroadcastsd 2760(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 90
               vbroadcastsd 2768(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 91
               vmovaps 5696(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 91
               vmovaps 5728(pA), ra1
            #endif
         #endif
         #if KB > 90
            vbroadcastsd 2776(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 91
               vbroadcastsd 2784(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 91
               vbroadcastsd 2792(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 91
               vbroadcastsd 2800(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 92
               vmovaps 5760(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 92
               vmovaps 5792(pA), rA1
            #endif
         #endif
         #if KB > 91
            vbroadcastsd 2808(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 92
               vbroadcastsd 2816(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 92
               vbroadcastsd 2824(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 92
               vbroadcastsd 2832(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 93
               vmovaps 5824(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 93
               vmovaps 5856(pA), ra1
            #endif
         #endif
         #if KB > 92
            vbroadcastsd 2840(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 93
               vbroadcastsd 2848(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 93
               vbroadcastsd 2856(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 93
               vbroadcastsd 2864(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 94
               vmovaps 5888(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 94
               vmovaps 5920(pA), rA1
            #endif
         #endif
         #if KB > 93
            vbroadcastsd 2872(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 94
               vbroadcastsd 2880(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 94
               vbroadcastsd 2888(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 94
               vbroadcastsd 2896(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 95
               vmovaps 5952(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 95
               vmovaps 5984(pA), ra1
            #endif
         #endif
         #if KB > 94
            vbroadcastsd 2904(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 95
               vbroadcastsd 2912(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 95
               vbroadcastsd 2920(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 95
               vbroadcastsd 2928(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 96
               vmovaps 6016(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 96
               vmovaps 6048(pA), rA1
            #endif
         #endif
         #if KB > 95
            vbroadcastsd 2936(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 96
               vbroadcastsd 2944(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 96
               vbroadcastsd 2952(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 96
               vbroadcastsd 2960(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 97
               vmovaps 6080(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 97
               vmovaps 6112(pA), ra1
            #endif
         #endif
         #if KB > 96
            vbroadcastsd 2968(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 97
               vbroadcastsd 2976(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 97
               vbroadcastsd 2984(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 97
               vbroadcastsd 2992(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 98
               vmovaps 6144(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 98
               vmovaps 6176(pA), rA1
            #endif
         #endif
         #if KB > 97
            vbroadcastsd 3000(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 98
               vbroadcastsd 3008(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 98
               vbroadcastsd 3016(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 98
               vbroadcastsd 3024(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 99
               vmovaps 6208(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 99
               vmovaps 6240(pA), ra1
            #endif
         #endif
         #if KB > 98
            vbroadcastsd 3032(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 99
               vbroadcastsd 3040(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 99
               vbroadcastsd 3048(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 99
               vbroadcastsd 3056(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 100
               vmovaps 6272(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 100
               vmovaps 6304(pA), rA1
            #endif
         #endif
         #if KB > 99
            vbroadcastsd 3064(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 100
               vbroadcastsd 3072(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 100
               vbroadcastsd 3080(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 100
               vbroadcastsd 3088(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 101
               vmovaps 6336(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 101
               vmovaps 6368(pA), ra1
            #endif
         #endif
         #if KB > 100
            vbroadcastsd 3096(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 101
               vbroadcastsd 3104(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 101
               vbroadcastsd 3112(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 101
               vbroadcastsd 3120(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 102
               vmovaps 6400(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 102
               vmovaps 6432(pA), rA1
            #endif
         #endif
         #if KB > 101
            vbroadcastsd 3128(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 102
               vbroadcastsd 3136(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 102
               vbroadcastsd 3144(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 102
               vbroadcastsd 3152(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 103
               vmovaps 6464(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 103
               vmovaps 6496(pA), ra1
            #endif
         #endif
         #if KB > 102
            vbroadcastsd 3160(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 103
               vbroadcastsd 3168(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 103
               vbroadcastsd 3176(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 103
               vbroadcastsd 3184(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 104
               vmovaps 6528(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 104
               vmovaps 6560(pA), rA1
            #endif
         #endif
         #if KB > 103
            vbroadcastsd 3192(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 104
               vbroadcastsd 3200(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 104
               vbroadcastsd 3208(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 104
               vbroadcastsd 3216(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 105
               vmovaps 6592(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 105
               vmovaps 6624(pA), ra1
            #endif
         #endif
         #if KB > 104
            vbroadcastsd 3224(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 105
               vbroadcastsd 3232(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 105
               vbroadcastsd 3240(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 105
               vbroadcastsd 3248(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 106
               vmovaps 6656(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 106
               vmovaps 6688(pA), rA1
            #endif
         #endif
         #if KB > 105
            vbroadcastsd 3256(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 106
               vbroadcastsd 3264(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 106
               vbroadcastsd 3272(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 106
               vbroadcastsd 3280(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 107
               vmovaps 6720(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 107
               vmovaps 6752(pA), ra1
            #endif
         #endif
         #if KB > 106
            vbroadcastsd 3288(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 107
               vbroadcastsd 3296(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 107
               vbroadcastsd 3304(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 107
               vbroadcastsd 3312(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 108
               vmovaps 6784(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 108
               vmovaps 6816(pA), rA1
            #endif
         #endif
         #if KB > 107
            vbroadcastsd 3320(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 108
               vbroadcastsd 3328(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 108
               vbroadcastsd 3336(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 108
               vbroadcastsd 3344(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 109
               vmovaps 6848(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 109
               vmovaps 6880(pA), ra1
            #endif
         #endif
         #if KB > 108
            vbroadcastsd 3352(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 109
               vbroadcastsd 3360(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 109
               vbroadcastsd 3368(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 109
               vbroadcastsd 3376(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 110
               vmovaps 6912(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 110
               vmovaps 6944(pA), rA1
            #endif
         #endif
         #if KB > 109
            vbroadcastsd 3384(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 110
               vbroadcastsd 3392(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 110
               vbroadcastsd 3400(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 110
               vbroadcastsd 3408(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 111
               vmovaps 6976(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 111
               vmovaps 7008(pA), ra1
            #endif
         #endif
         #if KB > 110
            vbroadcastsd 3416(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 111
               vbroadcastsd 3424(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 111
               vbroadcastsd 3432(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 111
               vbroadcastsd 3440(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 112
               vmovaps 7040(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 112
               vmovaps 7072(pA), rA1
            #endif
         #endif
         #if KB > 111
            vbroadcastsd 3448(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 112
               vbroadcastsd 3456(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 112
               vbroadcastsd 3464(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 112
               vbroadcastsd 3472(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 113
               vmovaps 7104(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 113
               vmovaps 7136(pA), ra1
            #endif
         #endif
         #if KB > 112
            vbroadcastsd 3480(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 113
               vbroadcastsd 3488(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 113
               vbroadcastsd 3496(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 113
               vbroadcastsd 3504(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 114
               vmovaps 7168(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 114
               vmovaps 7200(pA), rA1
            #endif
         #endif
         #if KB > 113
            vbroadcastsd 3512(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 114
               vbroadcastsd 3520(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 114
               vbroadcastsd 3528(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 114
               vbroadcastsd 3536(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 115
               vmovaps 7232(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 115
               vmovaps 7264(pA), ra1
            #endif
         #endif
         #if KB > 114
            vbroadcastsd 3544(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 115
               vbroadcastsd 3552(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 115
               vbroadcastsd 3560(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 115
               vbroadcastsd 3568(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 116
               vmovaps 7296(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 116
               vmovaps 7328(pA), rA1
            #endif
         #endif
         #if KB > 115
            vbroadcastsd 3576(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 116
               vbroadcastsd 3584(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 116
               vbroadcastsd 3592(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 116
               vbroadcastsd 3600(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 117
               vmovaps 7360(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 117
               vmovaps 7392(pA), ra1
            #endif
         #endif
         #if KB > 116
            vbroadcastsd 3608(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 117
               vbroadcastsd 3616(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 117
               vbroadcastsd 3624(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 117
               vbroadcastsd 3632(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 118
               vmovaps 7424(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 118
               vmovaps 7456(pA), rA1
            #endif
         #endif
         #if KB > 117
            vbroadcastsd 3640(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 118
               vbroadcastsd 3648(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 118
               vbroadcastsd 3656(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 118
               vbroadcastsd 3664(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 119
               vmovaps 7488(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 119
               vmovaps 7520(pA), ra1
            #endif
         #endif
         #if KB > 118
            vbroadcastsd 3672(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 119
               vbroadcastsd 3680(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 119
               vbroadcastsd 3688(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 119
               vbroadcastsd 3696(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 120
               vmovaps 7552(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 120
               vmovaps 7584(pA), rA1
            #endif
         #endif
         #if KB > 119
            vbroadcastsd 3704(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 120
               vbroadcastsd 3712(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 120
               vbroadcastsd 3720(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 120
               vbroadcastsd 3728(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 121
               vmovaps 7616(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 121
               vmovaps 7648(pA), ra1
            #endif
         #endif
         #if KB > 120
            vbroadcastsd 3736(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 121
               vbroadcastsd 3744(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 121
               vbroadcastsd 3752(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 121
               vbroadcastsd 3760(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 122
               vmovaps 7680(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 122
               vmovaps 7712(pA), rA1
            #endif
         #endif
         #if KB > 121
            vbroadcastsd 3768(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 122
               vbroadcastsd 3776(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 122
               vbroadcastsd 3784(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 122
               vbroadcastsd 3792(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 123
               vmovaps 7744(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 123
               vmovaps 7776(pA), ra1
            #endif
         #endif
         #if KB > 122
            vbroadcastsd 3800(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 123
               vbroadcastsd 3808(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 123
               vbroadcastsd 3816(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 123
               vbroadcastsd 3824(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 124
               vmovaps 7808(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 124
               vmovaps 7840(pA), rA1
            #endif
         #endif
         #if KB > 123
            vbroadcastsd 3832(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 124
               vbroadcastsd 3840(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 124
               vbroadcastsd 3848(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 124
               vbroadcastsd 3856(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 125
               vmovaps 7872(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 125
               vmovaps 7904(pA), ra1
            #endif
         #endif
         #if KB > 124
            vbroadcastsd 3864(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 125
               vbroadcastsd 3872(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 125
               vbroadcastsd 3880(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 125
               vbroadcastsd 3888(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 126
               vmovaps 7936(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 126
               vmovaps 7968(pA), rA1
            #endif
         #endif
         #if KB > 125
            vbroadcastsd 3896(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 126
               vbroadcastsd 3904(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 126
               vbroadcastsd 3912(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 126
               vbroadcastsd 3920(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 127
               vmovaps 8000(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 127
               vmovaps 8032(pA), ra1
            #endif
         #endif
         #if KB > 126
            vbroadcastsd 3928(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 127
               vbroadcastsd 3936(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 127
               vbroadcastsd 3944(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 127
               vbroadcastsd 3952(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 128
               vmovaps 8064(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 128
               vmovaps 8096(pA), rA1
            #endif
         #endif
         #if KB > 127
            vbroadcastsd 3960(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 128
               vbroadcastsd 3968(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 128
               vbroadcastsd 3976(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 128
               vbroadcastsd 3984(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 129
               vmovaps 8128(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 129
               vmovaps 8160(pA), ra1
            #endif
         #endif
         #if KB > 128
            vbroadcastsd 3992(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 129
               vbroadcastsd 4000(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 129
               vbroadcastsd 4008(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 129
               vbroadcastsd 4016(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 130
               vmovaps 8192(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 130
               vmovaps 8224(pA), rA1
            #endif
         #endif
         #if KB > 129
            vbroadcastsd 4024(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 130
               vbroadcastsd 4032(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 130
               vbroadcastsd 4040(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 130
               vbroadcastsd 4048(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 131
               vmovaps 8256(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 131
               vmovaps 8288(pA), ra1
            #endif
         #endif
         #if KB > 130
            vbroadcastsd 4056(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 131
               vbroadcastsd 4064(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 131
               vbroadcastsd 4072(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 131
               vbroadcastsd 4080(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 132
               vmovaps 8320(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 132
               vmovaps 8352(pA), rA1
            #endif
         #endif
         #if KB > 131
            vbroadcastsd 4088(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 132
               vbroadcastsd 4096(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 132
               vbroadcastsd 4104(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 132
               vbroadcastsd 4112(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 133
               vmovaps 8384(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 133
               vmovaps 8416(pA), ra1
            #endif
         #endif
         #if KB > 132
            vbroadcastsd 4120(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 133
               vbroadcastsd 4128(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 133
               vbroadcastsd 4136(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 133
               vbroadcastsd 4144(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 134
               vmovaps 8448(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 134
               vmovaps 8480(pA), rA1
            #endif
         #endif
         #if KB > 133
            vbroadcastsd 4152(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 134
               vbroadcastsd 4160(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 134
               vbroadcastsd 4168(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 134
               vbroadcastsd 4176(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 135
               vmovaps 8512(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 135
               vmovaps 8544(pA), ra1
            #endif
         #endif
         #if KB > 134
            vbroadcastsd 4184(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 135
               vbroadcastsd 4192(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 135
               vbroadcastsd 4200(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 135
               vbroadcastsd 4208(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 136
               vmovaps 8576(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 136
               vmovaps 8608(pA), rA1
            #endif
         #endif
         #if KB > 135
            vbroadcastsd 4216(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 136
               vbroadcastsd 4224(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 136
               vbroadcastsd 4232(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 136
               vbroadcastsd 4240(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 137
               vmovaps 8640(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 137
               vmovaps 8672(pA), ra1
            #endif
         #endif
         #if KB > 136
            vbroadcastsd 4248(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 137
               vbroadcastsd 4256(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 137
               vbroadcastsd 4264(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 137
               vbroadcastsd 4272(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 138
               vmovaps 8704(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 138
               vmovaps 8736(pA), rA1
            #endif
         #endif
         #if KB > 137
            vbroadcastsd 4280(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 138
               vbroadcastsd 4288(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 138
               vbroadcastsd 4296(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 138
               vbroadcastsd 4304(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 139
               vmovaps 8768(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 139
               vmovaps 8800(pA), ra1
            #endif
         #endif
         #if KB > 138
            vbroadcastsd 4312(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 139
               vbroadcastsd 4320(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 139
               vbroadcastsd 4328(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 139
               vbroadcastsd 4336(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 140
               vmovaps 8832(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 140
               vmovaps 8864(pA), rA1
            #endif
         #endif
         #if KB > 139
            vbroadcastsd 4344(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 140
               vbroadcastsd 4352(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 140
               vbroadcastsd 4360(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 140
               vbroadcastsd 4368(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 141
               vmovaps 8896(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 141
               vmovaps 8928(pA), ra1
            #endif
         #endif
         #if KB > 140
            vbroadcastsd 4376(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 141
               vbroadcastsd 4384(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 141
               vbroadcastsd 4392(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 141
               vbroadcastsd 4400(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 142
               vmovaps 8960(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 142
               vmovaps 8992(pA), rA1
            #endif
         #endif
         #if KB > 141
            vbroadcastsd 4408(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 142
               vbroadcastsd 4416(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 142
               vbroadcastsd 4424(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 142
               vbroadcastsd 4432(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 143
               vmovaps 9024(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 143
               vmovaps 9056(pA), ra1
            #endif
         #endif
         #if KB > 142
            vbroadcastsd 4440(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 143
               vbroadcastsd 4448(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 143
               vbroadcastsd 4456(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 143
               vbroadcastsd 4464(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 144
               vmovaps 9088(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 144
               vmovaps 9120(pA), rA1
            #endif
         #endif
         #if KB > 143
            vbroadcastsd 4472(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 144
               vbroadcastsd 4480(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 144
               vbroadcastsd 4488(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 144
               vbroadcastsd 4496(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 145
               vmovaps 9152(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 145
               vmovaps 9184(pA), ra1
            #endif
         #endif
         #if KB > 144
            vbroadcastsd 4504(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 145
               vbroadcastsd 4512(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 145
               vbroadcastsd 4520(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 145
               vbroadcastsd 4528(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 146
               vmovaps 9216(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 146
               vmovaps 9248(pA), rA1
            #endif
         #endif
         #if KB > 145
            vbroadcastsd 4536(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 146
               vbroadcastsd 4544(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 146
               vbroadcastsd 4552(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 146
               vbroadcastsd 4560(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 147
               vmovaps 9280(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 147
               vmovaps 9312(pA), ra1
            #endif
         #endif
         #if KB > 146
            vbroadcastsd 4568(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 147
               vbroadcastsd 4576(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 147
               vbroadcastsd 4584(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 147
               vbroadcastsd 4592(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 148
               vmovaps 9344(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 148
               vmovaps 9376(pA), rA1
            #endif
         #endif
         #if KB > 147
            vbroadcastsd 4600(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 148
               vbroadcastsd 4608(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 148
               vbroadcastsd 4616(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 148
               vbroadcastsd 4624(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 149
               vmovaps 9408(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 149
               vmovaps 9440(pA), ra1
            #endif
         #endif
         #if KB > 148
            vbroadcastsd 4632(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 149
               vbroadcastsd 4640(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 149
               vbroadcastsd 4648(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 149
               vbroadcastsd 4656(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 150
               vmovaps 9472(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 150
               vmovaps 9504(pA), rA1
            #endif
         #endif
         #if KB > 149
            vbroadcastsd 4664(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 150
               vbroadcastsd 4672(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 150
               vbroadcastsd 4680(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 150
               vbroadcastsd 4688(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 151
               vmovaps 9536(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 151
               vmovaps 9568(pA), ra1
            #endif
         #endif
         #if KB > 150
            vbroadcastsd 4696(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 151
               vbroadcastsd 4704(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 151
               vbroadcastsd 4712(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 151
               vbroadcastsd 4720(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 152
               vmovaps 9600(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 152
               vmovaps 9632(pA), rA1
            #endif
         #endif
         #if KB > 151
            vbroadcastsd 4728(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 152
               vbroadcastsd 4736(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 152
               vbroadcastsd 4744(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 152
               vbroadcastsd 4752(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 153
               vmovaps 9664(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 153
               vmovaps 9696(pA), ra1
            #endif
         #endif
         #if KB > 152
            vbroadcastsd 4760(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 153
               vbroadcastsd 4768(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 153
               vbroadcastsd 4776(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 153
               vbroadcastsd 4784(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 154
               vmovaps 9728(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 154
               vmovaps 9760(pA), rA1
            #endif
         #endif
         #if KB > 153
            vbroadcastsd 4792(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 154
               vbroadcastsd 4800(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 154
               vbroadcastsd 4808(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 154
               vbroadcastsd 4816(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 155
               vmovaps 9792(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 155
               vmovaps 9824(pA), ra1
            #endif
         #endif
         #if KB > 154
            vbroadcastsd 4824(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 155
               vbroadcastsd 4832(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 155
               vbroadcastsd 4840(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 155
               vbroadcastsd 4848(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 156
               vmovaps 9856(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 156
               vmovaps 9888(pA), rA1
            #endif
         #endif
         #if KB > 155
            vbroadcastsd 4856(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 156
               vbroadcastsd 4864(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 156
               vbroadcastsd 4872(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 156
               vbroadcastsd 4880(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 157
               vmovaps 9920(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 157
               vmovaps 9952(pA), ra1
            #endif
         #endif
         #if KB > 156
            vbroadcastsd 4888(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 157
               vbroadcastsd 4896(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 157
               vbroadcastsd 4904(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 157
               vbroadcastsd 4912(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 158
               vmovaps 9984(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 158
               vmovaps 10016(pA), rA1
            #endif
         #endif
         #if KB > 157
            vbroadcastsd 4920(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 158
               vbroadcastsd 4928(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 158
               vbroadcastsd 4936(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 158
               vbroadcastsd 4944(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 159
               vmovaps 10048(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 159
               vmovaps 10080(pA), ra1
            #endif
         #endif
         #if KB > 158
            vbroadcastsd 4952(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 159
               vbroadcastsd 4960(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 159
               vbroadcastsd 4968(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 159
               vbroadcastsd 4976(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 160
               vmovaps 10112(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 160
               vmovaps 10144(pA), rA1
            #endif
         #endif
         #if KB > 159
            vbroadcastsd 4984(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 160
               vbroadcastsd 4992(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 160
               vbroadcastsd 5000(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 160
               vbroadcastsd 5008(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 161
               vmovaps 10176(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 161
               vmovaps 10208(pA), ra1
            #endif
         #endif
         #if KB > 160
            vbroadcastsd 5016(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 161
               vbroadcastsd 5024(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 161
               vbroadcastsd 5032(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 161
               vbroadcastsd 5040(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 162
               vmovaps 10240(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 162
               vmovaps 10272(pA), rA1
            #endif
         #endif
         #if KB > 161
            vbroadcastsd 5048(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 162
               vbroadcastsd 5056(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 162
               vbroadcastsd 5064(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 162
               vbroadcastsd 5072(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 163
               vmovaps 10304(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 163
               vmovaps 10336(pA), ra1
            #endif
         #endif
         #if KB > 162
            vbroadcastsd 5080(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 163
               vbroadcastsd 5088(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 163
               vbroadcastsd 5096(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 163
               vbroadcastsd 5104(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 164
               vmovaps 10368(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 164
               vmovaps 10400(pA), rA1
            #endif
         #endif
         #if KB > 163
            vbroadcastsd 5112(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 164
               vbroadcastsd 5120(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 164
               vbroadcastsd 5128(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 164
               vbroadcastsd 5136(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 165
               vmovaps 10432(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 165
               vmovaps 10464(pA), ra1
            #endif
         #endif
         #if KB > 164
            vbroadcastsd 5144(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 165
               vbroadcastsd 5152(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 165
               vbroadcastsd 5160(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 165
               vbroadcastsd 5168(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 166
               vmovaps 10496(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 166
               vmovaps 10528(pA), rA1
            #endif
         #endif
         #if KB > 165
            vbroadcastsd 5176(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 166
               vbroadcastsd 5184(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 166
               vbroadcastsd 5192(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 166
               vbroadcastsd 5200(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 167
               vmovaps 10560(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 167
               vmovaps 10592(pA), ra1
            #endif
         #endif
         #if KB > 166
            vbroadcastsd 5208(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 167
               vbroadcastsd 5216(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 167
               vbroadcastsd 5224(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 167
               vbroadcastsd 5232(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 168
               vmovaps 10624(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 168
               vmovaps 10656(pA), rA1
            #endif
         #endif
         #if KB > 167
            vbroadcastsd 5240(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 168
               vbroadcastsd 5248(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 168
               vbroadcastsd 5256(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 168
               vbroadcastsd 5264(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 169
               vmovaps 10688(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 169
               vmovaps 10720(pA), ra1
            #endif
         #endif
         #if KB > 168
            vbroadcastsd 5272(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 169
               vbroadcastsd 5280(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 169
               vbroadcastsd 5288(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 169
               vbroadcastsd 5296(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 170
               vmovaps 10752(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 170
               vmovaps 10784(pA), rA1
            #endif
         #endif
         #if KB > 169
            vbroadcastsd 5304(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 170
               vbroadcastsd 5312(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 170
               vbroadcastsd 5320(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 170
               vbroadcastsd 5328(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 171
               vmovaps 10816(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 171
               vmovaps 10848(pA), ra1
            #endif
         #endif
         #if KB > 170
            vbroadcastsd 5336(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 171
               vbroadcastsd 5344(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 171
               vbroadcastsd 5352(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 171
               vbroadcastsd 5360(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 172
               vmovaps 10880(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 172
               vmovaps 10912(pA), rA1
            #endif
         #endif
         #if KB > 171
            vbroadcastsd 5368(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 172
               vbroadcastsd 5376(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 172
               vbroadcastsd 5384(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 172
               vbroadcastsd 5392(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 173
               vmovaps 10944(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 173
               vmovaps 10976(pA), ra1
            #endif
         #endif
         #if KB > 172
            vbroadcastsd 5400(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 173
               vbroadcastsd 5408(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 173
               vbroadcastsd 5416(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 173
               vbroadcastsd 5424(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 174
               vmovaps 11008(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 174
               vmovaps 11040(pA), rA1
            #endif
         #endif
         #if KB > 173
            vbroadcastsd 5432(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 174
               vbroadcastsd 5440(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 174
               vbroadcastsd 5448(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 174
               vbroadcastsd 5456(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 175
               vmovaps 11072(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 175
               vmovaps 11104(pA), ra1
            #endif
         #endif
         #if KB > 174
            vbroadcastsd 5464(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 175
               vbroadcastsd 5472(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 175
               vbroadcastsd 5480(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 175
               vbroadcastsd 5488(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 176
               vmovaps 11136(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 176
               vmovaps 11168(pA), rA1
            #endif
         #endif
         #if KB > 175
            vbroadcastsd 5496(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 176
               vbroadcastsd 5504(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 176
               vbroadcastsd 5512(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 176
               vbroadcastsd 5520(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 177
               vmovaps 11200(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 177
               vmovaps 11232(pA), ra1
            #endif
         #endif
         #if KB > 176
            vbroadcastsd 5528(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 177
               vbroadcastsd 5536(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 177
               vbroadcastsd 5544(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 177
               vbroadcastsd 5552(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 178
               vmovaps 11264(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 178
               vmovaps 11296(pA), rA1
            #endif
         #endif
         #if KB > 177
            vbroadcastsd 5560(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 178
               vbroadcastsd 5568(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 178
               vbroadcastsd 5576(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 178
               vbroadcastsd 5584(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 179
               vmovaps 11328(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 179
               vmovaps 11360(pA), ra1
            #endif
         #endif
         #if KB > 178
            vbroadcastsd 5592(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 179
               vbroadcastsd 5600(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 179
               vbroadcastsd 5608(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 179
               vbroadcastsd 5616(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 180
               vmovaps 11392(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 180
               vmovaps 11424(pA), rA1
            #endif
         #endif
         #if KB > 179
            vbroadcastsd 5624(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 180
               vbroadcastsd 5632(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 180
               vbroadcastsd 5640(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 180
               vbroadcastsd 5648(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 181
               vmovaps 11456(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 181
               vmovaps 11488(pA), ra1
            #endif
         #endif
         #if KB > 180
            vbroadcastsd 5656(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 181
               vbroadcastsd 5664(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 181
               vbroadcastsd 5672(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 181
               vbroadcastsd 5680(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 182
               vmovaps 11520(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 182
               vmovaps 11552(pA), rA1
            #endif
         #endif
         #if KB > 181
            vbroadcastsd 5688(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 182
               vbroadcastsd 5696(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 182
               vbroadcastsd 5704(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 182
               vbroadcastsd 5712(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 183
               vmovaps 11584(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 183
               vmovaps 11616(pA), ra1
            #endif
         #endif
         #if KB > 182
            vbroadcastsd 5720(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 183
               vbroadcastsd 5728(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 183
               vbroadcastsd 5736(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 183
               vbroadcastsd 5744(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 184
               vmovaps 11648(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 184
               vmovaps 11680(pA), rA1
            #endif
         #endif
         #if KB > 183
            vbroadcastsd 5752(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 184
               vbroadcastsd 5760(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 184
               vbroadcastsd 5768(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 184
               vbroadcastsd 5776(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 185
               vmovaps 11712(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 185
               vmovaps 11744(pA), ra1
            #endif
         #endif
         #if KB > 184
            vbroadcastsd 5784(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 185
               vbroadcastsd 5792(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 185
               vbroadcastsd 5800(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 185
               vbroadcastsd 5808(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 186
               vmovaps 11776(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 186
               vmovaps 11808(pA), rA1
            #endif
         #endif
         #if KB > 185
            vbroadcastsd 5816(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 186
               vbroadcastsd 5824(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 186
               vbroadcastsd 5832(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 186
               vbroadcastsd 5840(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 187
               vmovaps 11840(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 187
               vmovaps 11872(pA), ra1
            #endif
         #endif
         #if KB > 186
            vbroadcastsd 5848(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 187
               vbroadcastsd 5856(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 187
               vbroadcastsd 5864(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 187
               vbroadcastsd 5872(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 188
               vmovaps 11904(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 188
               vmovaps 11936(pA), rA1
            #endif
         #endif
         #if KB > 187
            vbroadcastsd 5880(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 188
               vbroadcastsd 5888(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 188
               vbroadcastsd 5896(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 188
               vbroadcastsd 5904(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 189
               vmovaps 11968(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 189
               vmovaps 12000(pA), ra1
            #endif
         #endif
         #if KB > 188
            vbroadcastsd 5912(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 189
               vbroadcastsd 5920(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 189
               vbroadcastsd 5928(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 189
               vbroadcastsd 5936(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 190
               vmovaps 12032(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 190
               vmovaps 12064(pA), rA1
            #endif
         #endif
         #if KB > 189
            vbroadcastsd 5944(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 190
               vbroadcastsd 5952(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 190
               vbroadcastsd 5960(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 190
               vbroadcastsd 5968(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 191
               vmovaps 12096(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 191
               vmovaps 12128(pA), ra1
            #endif
         #endif
         #if KB > 190
            vbroadcastsd 5976(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 191
               vbroadcastsd 5984(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 191
               vbroadcastsd 5992(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 191
               vbroadcastsd 6000(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 192
               vmovaps 12160(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 192
               vmovaps 12192(pA), rA1
            #endif
         #endif
         #if KB > 191
            vbroadcastsd 6008(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 192
               vbroadcastsd 6016(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 192
               vbroadcastsd 6024(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 192
               vbroadcastsd 6032(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 193
               vmovaps 12224(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 193
               vmovaps 12256(pA), ra1
            #endif
         #endif
         #if KB > 192
            vbroadcastsd 6040(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 193
               vbroadcastsd 6048(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 193
               vbroadcastsd 6056(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 193
               vbroadcastsd 6064(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 194
               vmovaps 12288(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 194
               vmovaps 12320(pA), rA1
            #endif
         #endif
         #if KB > 193
            vbroadcastsd 6072(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 194
               vbroadcastsd 6080(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 194
               vbroadcastsd 6088(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 194
               vbroadcastsd 6096(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 195
               vmovaps 12352(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 195
               vmovaps 12384(pA), ra1
            #endif
         #endif
         #if KB > 194
            vbroadcastsd 6104(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 195
               vbroadcastsd 6112(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 195
               vbroadcastsd 6120(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 195
               vbroadcastsd 6128(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 196
               vmovaps 12416(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 196
               vmovaps 12448(pA), rA1
            #endif
         #endif
         #if KB > 195
            vbroadcastsd 6136(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 196
               vbroadcastsd 6144(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 196
               vbroadcastsd 6152(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 196
               vbroadcastsd 6160(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 197
               vmovaps 12480(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 197
               vmovaps 12512(pA), ra1
            #endif
         #endif
         #if KB > 196
            vbroadcastsd 6168(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 197
               vbroadcastsd 6176(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 197
               vbroadcastsd 6184(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 197
               vbroadcastsd 6192(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 198
               vmovaps 12544(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 198
               vmovaps 12576(pA), rA1
            #endif
         #endif
         #if KB > 197
            vbroadcastsd 6200(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 198
               vbroadcastsd 6208(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 198
               vbroadcastsd 6216(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 198
               vbroadcastsd 6224(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 199
               vmovaps 12608(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 199
               vmovaps 12640(pA), ra1
            #endif
         #endif
         #if KB > 198
            vbroadcastsd 6232(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 199
               vbroadcastsd 6240(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 199
               vbroadcastsd 6248(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 199
               vbroadcastsd 6256(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 200
               vmovaps 12672(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 200
               vmovaps 12704(pA), rA1
            #endif
         #endif
         #if KB > 199
            vbroadcastsd 6264(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 200
               vbroadcastsd 6272(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 200
               vbroadcastsd 6280(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 200
               vbroadcastsd 6288(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 201
               vmovaps 12736(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 201
               vmovaps 12768(pA), ra1
            #endif
         #endif
         #if KB > 200
            vbroadcastsd 6296(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 201
               vbroadcastsd 6304(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 201
               vbroadcastsd 6312(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 201
               vbroadcastsd 6320(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 202
               vmovaps 12800(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 202
               vmovaps 12832(pA), rA1
            #endif
         #endif
         #if KB > 201
            vbroadcastsd 6328(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 202
               vbroadcastsd 6336(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 202
               vbroadcastsd 6344(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 202
               vbroadcastsd 6352(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 203
               vmovaps 12864(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 203
               vmovaps 12896(pA), ra1
            #endif
         #endif
         #if KB > 202
            vbroadcastsd 6360(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 203
               vbroadcastsd 6368(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 203
               vbroadcastsd 6376(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 203
               vbroadcastsd 6384(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 204
               vmovaps 12928(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 204
               vmovaps 12960(pA), rA1
            #endif
         #endif
         #if KB > 203
            vbroadcastsd 6392(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 204
               vbroadcastsd 6400(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 204
               vbroadcastsd 6408(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 204
               vbroadcastsd 6416(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 205
               vmovaps 12992(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 205
               vmovaps 13024(pA), ra1
            #endif
         #endif
         #if KB > 204
            vbroadcastsd 6424(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 205
               vbroadcastsd 6432(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 205
               vbroadcastsd 6440(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 205
               vbroadcastsd 6448(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 206
               vmovaps 13056(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 206
               vmovaps 13088(pA), rA1
            #endif
         #endif
         #if KB > 205
            vbroadcastsd 6456(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 206
               vbroadcastsd 6464(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 206
               vbroadcastsd 6472(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 206
               vbroadcastsd 6480(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 207
               vmovaps 13120(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 207
               vmovaps 13152(pA), ra1
            #endif
         #endif
         #if KB > 206
            vbroadcastsd 6488(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 207
               vbroadcastsd 6496(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 207
               vbroadcastsd 6504(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 207
               vbroadcastsd 6512(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 208
               vmovaps 13184(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 208
               vmovaps 13216(pA), rA1
            #endif
         #endif
         #if KB > 207
            vbroadcastsd 6520(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 208
               vbroadcastsd 6528(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 208
               vbroadcastsd 6536(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 208
               vbroadcastsd 6544(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 209
               vmovaps 13248(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 209
               vmovaps 13280(pA), ra1
            #endif
         #endif
         #if KB > 208
            vbroadcastsd 6552(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 209
               vbroadcastsd 6560(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 209
               vbroadcastsd 6568(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 209
               vbroadcastsd 6576(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 210
               vmovaps 13312(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 210
               vmovaps 13344(pA), rA1
            #endif
         #endif
         #if KB > 209
            vbroadcastsd 6584(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 210
               vbroadcastsd 6592(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 210
               vbroadcastsd 6600(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 210
               vbroadcastsd 6608(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 211
               vmovaps 13376(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 211
               vmovaps 13408(pA), ra1
            #endif
         #endif
         #if KB > 210
            vbroadcastsd 6616(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 211
               vbroadcastsd 6624(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 211
               vbroadcastsd 6632(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 211
               vbroadcastsd 6640(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 212
               vmovaps 13440(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 212
               vmovaps 13472(pA), rA1
            #endif
         #endif
         #if KB > 211
            vbroadcastsd 6648(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 212
               vbroadcastsd 6656(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 212
               vbroadcastsd 6664(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 212
               vbroadcastsd 6672(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 213
               vmovaps 13504(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 213
               vmovaps 13536(pA), ra1
            #endif
         #endif
         #if KB > 212
            vbroadcastsd 6680(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 213
               vbroadcastsd 6688(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 213
               vbroadcastsd 6696(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 213
               vbroadcastsd 6704(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 214
               vmovaps 13568(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 214
               vmovaps 13600(pA), rA1
            #endif
         #endif
         #if KB > 213
            vbroadcastsd 6712(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 214
               vbroadcastsd 6720(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 214
               vbroadcastsd 6728(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 214
               vbroadcastsd 6736(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 215
               vmovaps 13632(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 215
               vmovaps 13664(pA), ra1
            #endif
         #endif
         #if KB > 214
            vbroadcastsd 6744(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 215
               vbroadcastsd 6752(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 215
               vbroadcastsd 6760(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 215
               vbroadcastsd 6768(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 216
               vmovaps 13696(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 216
               vmovaps 13728(pA), rA1
            #endif
         #endif
         #if KB > 215
            vbroadcastsd 6776(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 216
               vbroadcastsd 6784(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 216
               vbroadcastsd 6792(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 216
               vbroadcastsd 6800(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 217
               vmovaps 13760(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 217
               vmovaps 13792(pA), ra1
            #endif
         #endif
         #if KB > 216
            vbroadcastsd 6808(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 217
               vbroadcastsd 6816(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 217
               vbroadcastsd 6824(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 217
               vbroadcastsd 6832(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 218
               vmovaps 13824(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 218
               vmovaps 13856(pA), rA1
            #endif
         #endif
         #if KB > 217
            vbroadcastsd 6840(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 218
               vbroadcastsd 6848(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 218
               vbroadcastsd 6856(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 218
               vbroadcastsd 6864(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 219
               vmovaps 13888(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 219
               vmovaps 13920(pA), ra1
            #endif
         #endif
         #if KB > 218
            vbroadcastsd 6872(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 219
               vbroadcastsd 6880(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 219
               vbroadcastsd 6888(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 219
               vbroadcastsd 6896(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 220
               vmovaps 13952(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 220
               vmovaps 13984(pA), rA1
            #endif
         #endif
         #if KB > 219
            vbroadcastsd 6904(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 220
               vbroadcastsd 6912(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 220
               vbroadcastsd 6920(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 220
               vbroadcastsd 6928(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 221
               vmovaps 14016(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 221
               vmovaps 14048(pA), ra1
            #endif
         #endif
         #if KB > 220
            vbroadcastsd 6936(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 221
               vbroadcastsd 6944(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 221
               vbroadcastsd 6952(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 221
               vbroadcastsd 6960(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 222
               vmovaps 14080(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 222
               vmovaps 14112(pA), rA1
            #endif
         #endif
         #if KB > 221
            vbroadcastsd 6968(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 222
               vbroadcastsd 6976(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 222
               vbroadcastsd 6984(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 222
               vbroadcastsd 6992(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 223
               vmovaps 14144(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 223
               vmovaps 14176(pA), ra1
            #endif
         #endif
         #if KB > 222
            vbroadcastsd 7000(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 223
               vbroadcastsd 7008(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 223
               vbroadcastsd 7016(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 223
               vbroadcastsd 7024(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 224
               vmovaps 14208(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 224
               vmovaps 14240(pA), rA1
            #endif
         #endif
         #if KB > 223
            vbroadcastsd 7032(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 224
               vbroadcastsd 7040(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 224
               vbroadcastsd 7048(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 224
               vbroadcastsd 7056(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 225
               vmovaps 14272(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 225
               vmovaps 14304(pA), ra1
            #endif
         #endif
         #if KB > 224
            vbroadcastsd 7064(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 225
               vbroadcastsd 7072(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 225
               vbroadcastsd 7080(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 225
               vbroadcastsd 7088(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 226
               vmovaps 14336(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 226
               vmovaps 14368(pA), rA1
            #endif
         #endif
         #if KB > 225
            vbroadcastsd 7096(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 226
               vbroadcastsd 7104(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 226
               vbroadcastsd 7112(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 226
               vbroadcastsd 7120(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 227
               vmovaps 14400(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 227
               vmovaps 14432(pA), ra1
            #endif
         #endif
         #if KB > 226
            vbroadcastsd 7128(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 227
               vbroadcastsd 7136(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 227
               vbroadcastsd 7144(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 227
               vbroadcastsd 7152(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 228
               vmovaps 14464(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 228
               vmovaps 14496(pA), rA1
            #endif
         #endif
         #if KB > 227
            vbroadcastsd 7160(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 228
               vbroadcastsd 7168(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 228
               vbroadcastsd 7176(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 228
               vbroadcastsd 7184(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 229
               vmovaps 14528(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 229
               vmovaps 14560(pA), ra1
            #endif
         #endif
         #if KB > 228
            vbroadcastsd 7192(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 229
               vbroadcastsd 7200(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 229
               vbroadcastsd 7208(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 229
               vbroadcastsd 7216(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 230
               vmovaps 14592(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 230
               vmovaps 14624(pA), rA1
            #endif
         #endif
         #if KB > 229
            vbroadcastsd 7224(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 230
               vbroadcastsd 7232(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 230
               vbroadcastsd 7240(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 230
               vbroadcastsd 7248(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 231
               vmovaps 14656(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 231
               vmovaps 14688(pA), ra1
            #endif
         #endif
         #if KB > 230
            vbroadcastsd 7256(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 231
               vbroadcastsd 7264(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 231
               vbroadcastsd 7272(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 231
               vbroadcastsd 7280(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 232
               vmovaps 14720(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 232
               vmovaps 14752(pA), rA1
            #endif
         #endif
         #if KB > 231
            vbroadcastsd 7288(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 232
               vbroadcastsd 7296(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 232
               vbroadcastsd 7304(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 232
               vbroadcastsd 7312(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 233
               vmovaps 14784(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 233
               vmovaps 14816(pA), ra1
            #endif
         #endif
         #if KB > 232
            vbroadcastsd 7320(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 233
               vbroadcastsd 7328(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 233
               vbroadcastsd 7336(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 233
               vbroadcastsd 7344(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 234
               vmovaps 14848(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 234
               vmovaps 14880(pA), rA1
            #endif
         #endif
         #if KB > 233
            vbroadcastsd 7352(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 234
               vbroadcastsd 7360(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 234
               vbroadcastsd 7368(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 234
               vbroadcastsd 7376(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 235
               vmovaps 14912(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 235
               vmovaps 14944(pA), ra1
            #endif
         #endif
         #if KB > 234
            vbroadcastsd 7384(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 235
               vbroadcastsd 7392(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 235
               vbroadcastsd 7400(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 235
               vbroadcastsd 7408(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 236
               vmovaps 14976(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 236
               vmovaps 15008(pA), rA1
            #endif
         #endif
         #if KB > 235
            vbroadcastsd 7416(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 236
               vbroadcastsd 7424(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 236
               vbroadcastsd 7432(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 236
               vbroadcastsd 7440(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 237
               vmovaps 15040(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 237
               vmovaps 15072(pA), ra1
            #endif
         #endif
         #if KB > 236
            vbroadcastsd 7448(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 237
               vbroadcastsd 7456(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 237
               vbroadcastsd 7464(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 237
               vbroadcastsd 7472(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 238
               vmovaps 15104(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 238
               vmovaps 15136(pA), rA1
            #endif
         #endif
         #if KB > 237
            vbroadcastsd 7480(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 238
               vbroadcastsd 7488(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 238
               vbroadcastsd 7496(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 238
               vbroadcastsd 7504(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 239
               vmovaps 15168(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 239
               vmovaps 15200(pA), ra1
            #endif
         #endif
         #if KB > 238
            vbroadcastsd 7512(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 239
               vbroadcastsd 7520(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 239
               vbroadcastsd 7528(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 239
               vbroadcastsd 7536(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 240
               vmovaps 15232(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 240
               vmovaps 15264(pA), rA1
            #endif
         #endif
         #if KB > 239
            vbroadcastsd 7544(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 240
               vbroadcastsd 7552(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 240
               vbroadcastsd 7560(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 240
               vbroadcastsd 7568(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 241
               vmovaps 15296(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 241
               vmovaps 15328(pA), ra1
            #endif
         #endif
         #if KB > 240
            vbroadcastsd 7576(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 241
               vbroadcastsd 7584(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 241
               vbroadcastsd 7592(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 241
               vbroadcastsd 7600(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 242
               vmovaps 15360(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 242
               vmovaps 15392(pA), rA1
            #endif
         #endif
         #if KB > 241
            vbroadcastsd 7608(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 242
               vbroadcastsd 7616(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 242
               vbroadcastsd 7624(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 242
               vbroadcastsd 7632(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 243
               vmovaps 15424(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 243
               vmovaps 15456(pA), ra1
            #endif
         #endif
         #if KB > 242
            vbroadcastsd 7640(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 243
               vbroadcastsd 7648(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 243
               vbroadcastsd 7656(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 243
               vbroadcastsd 7664(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 244
               vmovaps 15488(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 244
               vmovaps 15520(pA), rA1
            #endif
         #endif
         #if KB > 243
            vbroadcastsd 7672(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 244
               vbroadcastsd 7680(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 244
               vbroadcastsd 7688(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 244
               vbroadcastsd 7696(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 245
               vmovaps 15552(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 245
               vmovaps 15584(pA), ra1
            #endif
         #endif
         #if KB > 244
            vbroadcastsd 7704(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 245
               vbroadcastsd 7712(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 245
               vbroadcastsd 7720(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 245
               vbroadcastsd 7728(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 246
               vmovaps 15616(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 246
               vmovaps 15648(pA), rA1
            #endif
         #endif
         #if KB > 245
            vbroadcastsd 7736(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 246
               vbroadcastsd 7744(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 246
               vbroadcastsd 7752(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 246
               vbroadcastsd 7760(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 247
               vmovaps 15680(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 247
               vmovaps 15712(pA), ra1
            #endif
         #endif
         #if KB > 246
            vbroadcastsd 7768(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 247
               vbroadcastsd 7776(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 247
               vbroadcastsd 7784(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 247
               vbroadcastsd 7792(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 248
               vmovaps 15744(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 248
               vmovaps 15776(pA), rA1
            #endif
         #endif
         #if KB > 247
            vbroadcastsd 7800(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 248
               vbroadcastsd 7808(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 248
               vbroadcastsd 7816(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 248
               vbroadcastsd 7824(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 249
               vmovaps 15808(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 249
               vmovaps 15840(pA), ra1
            #endif
         #endif
         #if KB > 248
            vbroadcastsd 7832(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 249
               vbroadcastsd 7840(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 249
               vbroadcastsd 7848(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 249
               vbroadcastsd 7856(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 250
               vmovaps 15872(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 250
               vmovaps 15904(pA), rA1
            #endif
         #endif
         #if KB > 249
            vbroadcastsd 7864(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 250
               vbroadcastsd 7872(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 250
               vbroadcastsd 7880(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 250
               vbroadcastsd 7888(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 251
               vmovaps 15936(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 251
               vmovaps 15968(pA), ra1
            #endif
         #endif
         #if KB > 250
            vbroadcastsd 7896(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 251
               vbroadcastsd 7904(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 251
               vbroadcastsd 7912(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 251
               vbroadcastsd 7920(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 252
               vmovaps 16000(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 252
               vmovaps 16032(pA), rA1
            #endif
         #endif
         #if KB > 251
            vbroadcastsd 7928(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 252
               vbroadcastsd 7936(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 252
               vbroadcastsd 7944(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 252
               vbroadcastsd 7952(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 253
               vmovaps 16064(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 253
               vmovaps 16096(pA), ra1
            #endif
         #endif
         #if KB > 252
            vbroadcastsd 7960(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 253
               vbroadcastsd 7968(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 253
               vbroadcastsd 7976(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 253
               vbroadcastsd 7984(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 254
               vmovaps 16128(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 254
               vmovaps 16160(pA), rA1
            #endif
         #endif
         #if KB > 253
            vbroadcastsd 7992(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 254
               vbroadcastsd 8000(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 254
               vbroadcastsd 8008(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 254
               vbroadcastsd 8016(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 255
               vmovaps 16192(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 255
               vmovaps 16224(pA), ra1
            #endif
         #endif
         #if KB > 254
            vbroadcastsd 8024(pB), rB3
            FMAC rA0, rB0, rC00
            FMAC rA1, rB0, rC10

            #if KB > 255
               vbroadcastsd 8032(pB), rB0
            #endif
            FMAC rA0, rB1, rC01
            FMAC rA1, rB1, rC11

            #if KB > 255
               vbroadcastsd 8040(pB), rB1
            #endif
            FMAC rA0, rB2, rC02
            FMAC rA1, rB2, rC12
            #if KB > 255
               vbroadcastsd 8048(pB), rB2
            #endif

            FMAC rA0, rB3, rC03
            #if KB > 256
               vmovaps 16256(pA), rA0
            #endif
            FMAC rA1, rB3, rC13
            #if KB > 256
               vmovaps 16288(pA), rA1
            #endif
         #endif
         #if KB > 255
            vbroadcastsd 8056(pB), rB3
            FMAC ra0, rB0, rC00
            FMAC ra1, rB0, rC10

            #if KB > 256
               vbroadcastsd 8064(pB), rB0
            #endif
            FMAC ra0, rB1, rC01
            FMAC ra1, rB1, rC11

            #if KB > 256
               vbroadcastsd 8072(pB), rB1
            #endif
            FMAC ra0, rB2, rC02
            FMAC ra1, rB2, rC12
            #if KB > 256
               vbroadcastsd 8080(pB), rB2
            #endif

            FMAC ra0, rB3, rC03
            #if KB > 257
               vmovaps 16320(pA), ra0
            #endif
            FMAC ra1, rB3, rC13
            #if KB > 257
               vmovaps 16352(pA), ra1
            #endif
         #endif

         #ifndef BETA0
            VCOP -128(pC), rC00, rC00
         #endif
         vmovapd rC00, -128(pC)
         #ifndef BETA0
            VCOP -96(pC), rC10, rC10
         #endif
         vmovapd rC10, -96(pC)
         #ifndef BETA0
            VCOP -64(pC), rC01, rC01
         #endif
         vmovapd rC01, -64(pC)
         #ifndef BETA0
            VCOP -32(pC), rC11, rC11
         #endif
         vmovapd rC11, -32(pC)
         #ifndef BETA0
            VCOP (pC), rC02, rC02
         #endif
         vmovapd rC02, (pC)
         #ifndef BETA0
            VCOP 32(pC), rC12, rC12
         #endif
         vmovapd rC12, 32(pC)
         #ifndef BETA0
            VCOP 64(pC), rC03, rC03
         #endif
         vmovapd rC03, 64(pC)
         #ifndef BETA0
            VCOP 96(pC), rC13, rC13
         #endif
         vmovapd rC13, 96(pC)

         add $KB*4*8, pB
         add $256, pC
         sub $1, nnu
      jnz NLOOP
      mov nnu0, nnu
      add incAm, pA0
      mov pA0, pA
      mov pB0, pB
      sub $1, nmu
   jnz MLOOP
 DONE:
   movq -8(%rsp), %rbp
   movq -16(%rsp), %rbx
   movq -24(%rsp), %r12
   movq -32(%rsp), %r13
   movq -40(%rsp), %r14
   movq -48(%rsp), %r15
   ret
#if 0
.global findSize
findSize:
mov $SS1-SS0, %rax
ret
SS0:
SS1:
#endif
