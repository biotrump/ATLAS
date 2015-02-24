#include "atlas_asm.h"

#define rB0     %ymm0
#define rA0     %ymm1
#define rA1     %ymm2
#define rA2     %ymm3
#define rC00    %ymm4
#define rC10    %ymm5
#define rC20    %ymm6
#define rC01    %ymm7
#define rC11    %ymm8
#define rC21    %ymm9
#define rC02    %ymm10
#define rC12    %ymm11
#define rC22    %ymm12
#define rC03    %ymm13
#define rC13    %ymm14
#define rC23    %ymm15
#define rC00x   %xmm4
#define rC10x   %xmm5
#define rC20x   %xmm6
#define rC01x   %xmm7
#define rC11x   %xmm8
#define rC21x   %xmm9
#define rC02x   %xmm10
#define rC12x   %xmm11
#define rC22x   %xmm12
#define rC03x   %xmm13
#define rC13x   %xmm14
#define rC23x   %xmm15
/*
 * Prioritize original registers for inner-loop operations, but inc regs
 * can be anything w/o changing opcode size, so use new regs for those
 */
#define KK      %rdx  /* API reg */
#define pA      %rcx  /* API reg */
#define pB      %rax  /* comes in as r9 */
#define r192    %r9   /* set after mov r9 to pC () */
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
#define prefA(m_) prefetcht0 m_
#define prefB(m_) prefetcht0 m_
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
   cmp pfA, pB
   CMOVE pfB, pfA
   CMOVEq 24(%rsp), pfB
   mov KK, KK0
   sub $-128, pA
   sub $-128, pB
   sub $-128, pC
   sub $-128, pfA
   sub $-128, pfB
   mov $192, r192
   mov pA, pA0
   mov pB, pB0
   mov $12*3*8, incPF  /* incPF = MU*NU*sizeof = 288 */
/*
 * incAm = 12*sizeof*K = 12*8*K = 3*4*8*K = 3*K*32
 */
   lea (KK, KK,2), incAm   /* incAm = 3*K */
   shl $5, incAm           /* incAm = 32*3*K */
   vxorps rA2, rA2, rA2

   ALIGN16
   MLOOP:
      NLOOP:
/*
 *       First iteration peeled to handle init of rC
 */
         vbroadcastsd      -128(pB), rB0
         vmovaps -128(pA), rA0
         vmulpd rA0, rB0, rC00
         vmovaps -96(pA), rA1
         vmulpd rA1, rB0, rC10
         vmovaps -64(pA), rA2
         vmulpd rA2, rB0, rC20

         vbroadcastsd      -120(pB), rB0
         vmulpd rA0, rB0, rC01
            prefA(-128(pfA))
            prefA((pfA))
         vmulpd rA1, rB0, rC11
         vmulpd rA2, rB0, rC21

         vbroadcastsd      -112(pB), rB0
         vmulpd rA0, rB0, rC02
         vmulpd rA1, rB0, rC12
         vmulpd rA2, rB0, rC22

         vbroadcastsd      -104(pB), rB0
         vmulpd rA0, rB0, rC03
            vmovapd -32(pA), rA0
            prefC(-128(pC))
         vmulpd rA1, rB0, rC13
            vmovapd (pA), rA1
            prefC((pC))
         vmulpd rA2, rB0, rC23

/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 1
            vbroadcastsd -96(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 32(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd -88(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd -80(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd -72(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 2
               vmovaps 64(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 2
               vmovaps 96(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 2
            vbroadcastsd -64(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 128(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd -56(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd -48(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd -40(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 3
               vmovaps 160(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 3
               vmovaps 192(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 3
            vbroadcastsd -32(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 224(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd -24(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd -16(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd -8(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 4
               vmovaps 256(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 4
               vmovaps 288(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 4
            vbroadcastsd 0(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 320(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 8(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 16(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 24(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 5
               vmovaps 352(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 5
               vmovaps 384(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 5
            vbroadcastsd 32(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 416(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 40(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 48(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 56(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 6
               vmovaps 448(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 6
               vmovaps 480(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 6
            vbroadcastsd 64(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 512(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 72(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 80(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 88(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 7
               vmovaps 544(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 7
               vmovaps 576(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 7
            vbroadcastsd 96(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 608(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 104(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 112(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 120(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 8
               vmovaps 640(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 8
               vmovaps 672(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 8
            vbroadcastsd 128(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 704(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 136(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 144(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 152(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 9
               vmovaps 736(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 9
               vmovaps 768(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 9
            vbroadcastsd 160(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 800(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 168(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 176(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 184(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 10
               vmovaps 832(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 10
               vmovaps 864(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 10
            vbroadcastsd 192(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 896(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 200(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 208(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 216(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 11
               vmovaps 928(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 11
               vmovaps 960(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 11
            vbroadcastsd 224(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 992(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 232(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 240(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 248(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 12
               vmovaps 1024(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 12
               vmovaps 1056(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 12
            vbroadcastsd 256(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1088(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 264(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 272(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 280(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 13
               vmovaps 1120(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 13
               vmovaps 1152(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 13
            vbroadcastsd 288(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1184(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 296(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 304(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 312(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 14
               vmovaps 1216(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 14
               vmovaps 1248(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 14
            vbroadcastsd 320(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1280(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 328(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 336(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 344(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 15
               vmovaps 1312(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 15
               vmovaps 1344(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 15
            vbroadcastsd 352(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1376(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 360(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 368(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 376(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 16
               vmovaps 1408(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 16
               vmovaps 1440(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 16
            vbroadcastsd 384(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1472(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 392(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 400(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 408(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 17
               vmovaps 1504(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 17
               vmovaps 1536(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 17
            vbroadcastsd 416(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1568(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 424(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 432(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 440(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 18
               vmovaps 1600(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 18
               vmovaps 1632(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 18
            vbroadcastsd 448(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1664(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 456(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 464(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 472(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 19
               vmovaps 1696(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 19
               vmovaps 1728(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 19
            vbroadcastsd 480(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1760(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 488(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 496(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 504(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 20
               vmovaps 1792(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 20
               vmovaps 1824(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 20
            vbroadcastsd 512(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1856(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 520(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 528(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 536(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 21
               vmovaps 1888(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 21
               vmovaps 1920(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 21
            vbroadcastsd 544(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1952(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 552(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 560(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 568(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 22
               vmovaps 1984(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 22
               vmovaps 2016(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 22
            vbroadcastsd 576(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2048(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 584(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 592(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 600(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 23
               vmovaps 2080(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 23
               vmovaps 2112(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 23
            vbroadcastsd 608(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2144(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 616(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 624(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 632(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 24
               vmovaps 2176(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 24
               vmovaps 2208(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 24
            vbroadcastsd 640(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2240(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 648(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 656(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 664(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 25
               vmovaps 2272(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 25
               vmovaps 2304(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 25
            vbroadcastsd 672(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2336(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 680(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 688(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 696(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 26
               vmovaps 2368(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 26
               vmovaps 2400(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 26
            vbroadcastsd 704(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2432(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 712(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 720(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 728(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 27
               vmovaps 2464(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 27
               vmovaps 2496(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 27
            vbroadcastsd 736(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2528(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 744(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 752(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 760(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 28
               vmovaps 2560(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 28
               vmovaps 2592(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 28
            vbroadcastsd 768(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2624(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 776(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 784(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 792(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 29
               vmovaps 2656(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 29
               vmovaps 2688(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 29
            vbroadcastsd 800(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2720(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 808(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 816(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 824(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 30
               vmovaps 2752(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 30
               vmovaps 2784(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 30
            vbroadcastsd 832(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2816(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 840(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 848(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 856(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 31
               vmovaps 2848(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 31
               vmovaps 2880(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 31
            vbroadcastsd 864(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2912(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 872(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 880(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 888(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 32
               vmovaps 2944(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 32
               vmovaps 2976(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 32
            vbroadcastsd 896(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3008(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 904(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 912(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 920(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 33
               vmovaps 3040(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 33
               vmovaps 3072(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 33
            vbroadcastsd 928(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3104(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 936(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 944(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 952(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 34
               vmovaps 3136(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 34
               vmovaps 3168(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 34
            vbroadcastsd 960(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3200(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 968(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 976(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 984(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 35
               vmovaps 3232(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 35
               vmovaps 3264(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 35
            vbroadcastsd 992(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3296(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1000(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1008(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1016(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 36
               vmovaps 3328(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 36
               vmovaps 3360(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 36
            vbroadcastsd 1024(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3392(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1032(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1040(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1048(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 37
               vmovaps 3424(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 37
               vmovaps 3456(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 37
            vbroadcastsd 1056(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3488(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1064(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1072(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1080(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 38
               vmovaps 3520(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 38
               vmovaps 3552(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 38
            vbroadcastsd 1088(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3584(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1096(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1104(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1112(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 39
               vmovaps 3616(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 39
               vmovaps 3648(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 39
            vbroadcastsd 1120(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3680(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1128(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1136(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1144(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 40
               vmovaps 3712(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 40
               vmovaps 3744(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 40
            vbroadcastsd 1152(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3776(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1160(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1168(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1176(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 41
               vmovaps 3808(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 41
               vmovaps 3840(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 41
            vbroadcastsd 1184(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3872(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1192(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1200(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1208(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 42
               vmovaps 3904(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 42
               vmovaps 3936(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 42
            vbroadcastsd 1216(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3968(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1224(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1232(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1240(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 43
               vmovaps 4000(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 43
               vmovaps 4032(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 43
            vbroadcastsd 1248(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4064(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1256(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1264(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1272(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 44
               vmovaps 4096(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 44
               vmovaps 4128(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 44
            vbroadcastsd 1280(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4160(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1288(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1296(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1304(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 45
               vmovaps 4192(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 45
               vmovaps 4224(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 45
            vbroadcastsd 1312(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4256(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1320(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1328(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1336(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 46
               vmovaps 4288(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 46
               vmovaps 4320(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 46
            vbroadcastsd 1344(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4352(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1352(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1360(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1368(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 47
               vmovaps 4384(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 47
               vmovaps 4416(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 47
            vbroadcastsd 1376(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4448(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1384(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1392(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1400(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 48
               vmovaps 4480(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 48
               vmovaps 4512(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 48
            vbroadcastsd 1408(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4544(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1416(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1424(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1432(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 49
               vmovaps 4576(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 49
               vmovaps 4608(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 49
            vbroadcastsd 1440(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4640(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1448(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1456(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1464(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 50
               vmovaps 4672(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 50
               vmovaps 4704(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 50
            vbroadcastsd 1472(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4736(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1480(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1488(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1496(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 51
               vmovaps 4768(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 51
               vmovaps 4800(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 51
            vbroadcastsd 1504(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4832(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1512(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1520(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1528(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 52
               vmovaps 4864(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 52
               vmovaps 4896(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 52
            vbroadcastsd 1536(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4928(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1544(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1552(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1560(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 53
               vmovaps 4960(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 53
               vmovaps 4992(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 53
            vbroadcastsd 1568(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5024(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1576(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1584(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1592(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 54
               vmovaps 5056(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 54
               vmovaps 5088(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 54
            vbroadcastsd 1600(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5120(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1608(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1616(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1624(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 55
               vmovaps 5152(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 55
               vmovaps 5184(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 55
            vbroadcastsd 1632(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5216(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1640(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1648(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1656(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 56
               vmovaps 5248(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 56
               vmovaps 5280(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 56
            vbroadcastsd 1664(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5312(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1672(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1680(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1688(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 57
               vmovaps 5344(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 57
               vmovaps 5376(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 57
            vbroadcastsd 1696(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5408(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1704(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1712(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1720(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 58
               vmovaps 5440(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 58
               vmovaps 5472(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 58
            vbroadcastsd 1728(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5504(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1736(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1744(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1752(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 59
               vmovaps 5536(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 59
               vmovaps 5568(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 59
            vbroadcastsd 1760(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5600(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1768(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1776(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1784(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 60
               vmovaps 5632(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 60
               vmovaps 5664(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 60
            vbroadcastsd 1792(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5696(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1800(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1808(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1816(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 61
               vmovaps 5728(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 61
               vmovaps 5760(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 61
            vbroadcastsd 1824(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5792(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1832(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1840(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1848(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 62
               vmovaps 5824(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 62
               vmovaps 5856(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 62
            vbroadcastsd 1856(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5888(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1864(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1872(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1880(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 63
               vmovaps 5920(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 63
               vmovaps 5952(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 63
            vbroadcastsd 1888(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5984(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1896(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1904(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1912(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 64
               vmovaps 6016(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 64
               vmovaps 6048(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 64
            vbroadcastsd 1920(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6080(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1928(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1936(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1944(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 65
               vmovaps 6112(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 65
               vmovaps 6144(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 65
            vbroadcastsd 1952(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6176(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1960(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 1968(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 1976(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 66
               vmovaps 6208(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 66
               vmovaps 6240(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 66
            vbroadcastsd 1984(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6272(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 1992(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2000(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2008(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 67
               vmovaps 6304(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 67
               vmovaps 6336(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 67
            vbroadcastsd 2016(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6368(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2024(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2032(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2040(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 68
               vmovaps 6400(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 68
               vmovaps 6432(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 68
            vbroadcastsd 2048(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6464(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2056(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2064(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2072(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 69
               vmovaps 6496(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 69
               vmovaps 6528(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 69
            vbroadcastsd 2080(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6560(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2088(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2096(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2104(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 70
               vmovaps 6592(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 70
               vmovaps 6624(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 70
            vbroadcastsd 2112(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6656(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2120(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2128(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2136(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 71
               vmovaps 6688(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 71
               vmovaps 6720(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 71
            vbroadcastsd 2144(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6752(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2152(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2160(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2168(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 72
               vmovaps 6784(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 72
               vmovaps 6816(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 72
            vbroadcastsd 2176(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6848(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2184(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2192(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2200(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 73
               vmovaps 6880(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 73
               vmovaps 6912(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 73
            vbroadcastsd 2208(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6944(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2216(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2224(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2232(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 74
               vmovaps 6976(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 74
               vmovaps 7008(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 74
            vbroadcastsd 2240(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7040(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2248(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2256(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2264(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 75
               vmovaps 7072(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 75
               vmovaps 7104(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 75
            vbroadcastsd 2272(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7136(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2280(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2288(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2296(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 76
               vmovaps 7168(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 76
               vmovaps 7200(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 76
            vbroadcastsd 2304(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7232(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2312(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2320(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2328(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 77
               vmovaps 7264(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 77
               vmovaps 7296(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 77
            vbroadcastsd 2336(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7328(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2344(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2352(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2360(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 78
               vmovaps 7360(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 78
               vmovaps 7392(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 78
            vbroadcastsd 2368(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7424(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2376(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2384(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2392(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 79
               vmovaps 7456(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 79
               vmovaps 7488(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 79
            vbroadcastsd 2400(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7520(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2408(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2416(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2424(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 80
               vmovaps 7552(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 80
               vmovaps 7584(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 80
            vbroadcastsd 2432(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7616(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2440(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2448(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2456(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 81
               vmovaps 7648(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 81
               vmovaps 7680(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 81
            vbroadcastsd 2464(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7712(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2472(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2480(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2488(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 82
               vmovaps 7744(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 82
               vmovaps 7776(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 82
            vbroadcastsd 2496(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7808(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2504(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2512(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2520(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 83
               vmovaps 7840(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 83
               vmovaps 7872(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 83
            vbroadcastsd 2528(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7904(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2536(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2544(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2552(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 84
               vmovaps 7936(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 84
               vmovaps 7968(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 84
            vbroadcastsd 2560(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8000(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2568(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2576(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2584(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 85
               vmovaps 8032(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 85
               vmovaps 8064(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 85
            vbroadcastsd 2592(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8096(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2600(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2608(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2616(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 86
               vmovaps 8128(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 86
               vmovaps 8160(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 86
            vbroadcastsd 2624(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8192(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2632(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2640(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2648(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 87
               vmovaps 8224(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 87
               vmovaps 8256(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 87
            vbroadcastsd 2656(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8288(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2664(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2672(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2680(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 88
               vmovaps 8320(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 88
               vmovaps 8352(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 88
            vbroadcastsd 2688(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8384(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2696(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2704(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2712(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 89
               vmovaps 8416(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 89
               vmovaps 8448(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 89
            vbroadcastsd 2720(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8480(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2728(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2736(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2744(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 90
               vmovaps 8512(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 90
               vmovaps 8544(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 90
            vbroadcastsd 2752(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8576(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2760(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2768(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2776(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 91
               vmovaps 8608(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 91
               vmovaps 8640(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 91
            vbroadcastsd 2784(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8672(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2792(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2800(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2808(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 92
               vmovaps 8704(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 92
               vmovaps 8736(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 92
            vbroadcastsd 2816(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8768(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2824(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2832(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2840(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 93
               vmovaps 8800(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 93
               vmovaps 8832(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 93
            vbroadcastsd 2848(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8864(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2856(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2864(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2872(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 94
               vmovaps 8896(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 94
               vmovaps 8928(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 94
            vbroadcastsd 2880(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8960(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2888(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2896(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2904(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 95
               vmovaps 8992(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 95
               vmovaps 9024(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 95
            vbroadcastsd 2912(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9056(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2920(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2928(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2936(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 96
               vmovaps 9088(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 96
               vmovaps 9120(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 96
            vbroadcastsd 2944(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9152(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2952(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2960(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 2968(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 97
               vmovaps 9184(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 97
               vmovaps 9216(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 97
            vbroadcastsd 2976(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9248(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 2984(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 2992(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3000(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 98
               vmovaps 9280(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 98
               vmovaps 9312(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 98
            vbroadcastsd 3008(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9344(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3016(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3024(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3032(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 99
               vmovaps 9376(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 99
               vmovaps 9408(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 99
            vbroadcastsd 3040(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9440(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3048(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3056(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3064(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 100
               vmovaps 9472(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 100
               vmovaps 9504(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 100
            vbroadcastsd 3072(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9536(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3080(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3088(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3096(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 101
               vmovaps 9568(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 101
               vmovaps 9600(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 101
            vbroadcastsd 3104(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9632(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3112(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3120(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3128(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 102
               vmovaps 9664(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 102
               vmovaps 9696(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 102
            vbroadcastsd 3136(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9728(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3144(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3152(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3160(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 103
               vmovaps 9760(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 103
               vmovaps 9792(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 103
            vbroadcastsd 3168(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9824(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3176(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3184(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3192(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 104
               vmovaps 9856(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 104
               vmovaps 9888(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 104
            vbroadcastsd 3200(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9920(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3208(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3216(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3224(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 105
               vmovaps 9952(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 105
               vmovaps 9984(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 105
            vbroadcastsd 3232(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10016(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3240(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3248(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3256(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 106
               vmovaps 10048(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 106
               vmovaps 10080(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 106
            vbroadcastsd 3264(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10112(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3272(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3280(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3288(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 107
               vmovaps 10144(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 107
               vmovaps 10176(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 107
            vbroadcastsd 3296(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10208(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3304(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3312(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3320(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 108
               vmovaps 10240(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 108
               vmovaps 10272(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 108
            vbroadcastsd 3328(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10304(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3336(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3344(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3352(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 109
               vmovaps 10336(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 109
               vmovaps 10368(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 109
            vbroadcastsd 3360(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10400(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3368(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3376(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3384(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 110
               vmovaps 10432(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 110
               vmovaps 10464(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 110
            vbroadcastsd 3392(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10496(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3400(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3408(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3416(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 111
               vmovaps 10528(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 111
               vmovaps 10560(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 111
            vbroadcastsd 3424(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10592(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3432(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3440(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3448(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 112
               vmovaps 10624(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 112
               vmovaps 10656(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 112
            vbroadcastsd 3456(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10688(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3464(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3472(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3480(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 113
               vmovaps 10720(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 113
               vmovaps 10752(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 113
            vbroadcastsd 3488(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10784(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3496(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3504(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3512(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 114
               vmovaps 10816(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 114
               vmovaps 10848(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 114
            vbroadcastsd 3520(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10880(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3528(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3536(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3544(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 115
               vmovaps 10912(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 115
               vmovaps 10944(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 115
            vbroadcastsd 3552(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10976(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3560(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3568(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3576(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 116
               vmovaps 11008(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 116
               vmovaps 11040(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 116
            vbroadcastsd 3584(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11072(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3592(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3600(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3608(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 117
               vmovaps 11104(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 117
               vmovaps 11136(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 117
            vbroadcastsd 3616(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11168(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3624(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3632(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3640(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 118
               vmovaps 11200(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 118
               vmovaps 11232(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 118
            vbroadcastsd 3648(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11264(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3656(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3664(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3672(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 119
               vmovaps 11296(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 119
               vmovaps 11328(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 119
            vbroadcastsd 3680(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11360(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3688(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3696(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3704(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 120
               vmovaps 11392(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 120
               vmovaps 11424(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 120
            vbroadcastsd 3712(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11456(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3720(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3728(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3736(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 121
               vmovaps 11488(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 121
               vmovaps 11520(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 121
            vbroadcastsd 3744(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11552(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3752(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3760(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3768(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 122
               vmovaps 11584(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 122
               vmovaps 11616(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 122
            vbroadcastsd 3776(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11648(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3784(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3792(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3800(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 123
               vmovaps 11680(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 123
               vmovaps 11712(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 123
            vbroadcastsd 3808(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11744(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3816(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3824(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3832(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 124
               vmovaps 11776(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 124
               vmovaps 11808(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 124
            vbroadcastsd 3840(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11840(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3848(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3856(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3864(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 125
               vmovaps 11872(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 125
               vmovaps 11904(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 125
            vbroadcastsd 3872(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11936(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3880(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3888(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3896(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 126
               vmovaps 11968(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 126
               vmovaps 12000(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 126
            vbroadcastsd 3904(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12032(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3912(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3920(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3928(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 127
               vmovaps 12064(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 127
               vmovaps 12096(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 127
            vbroadcastsd 3936(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12128(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3944(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3952(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3960(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 128
               vmovaps 12160(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 128
               vmovaps 12192(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 128
            vbroadcastsd 3968(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12224(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 3976(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 3984(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 3992(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 129
               vmovaps 12256(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 129
               vmovaps 12288(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 129
            vbroadcastsd 4000(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12320(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4008(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4016(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4024(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 130
               vmovaps 12352(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 130
               vmovaps 12384(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 130
            vbroadcastsd 4032(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12416(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4040(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4048(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4056(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 131
               vmovaps 12448(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 131
               vmovaps 12480(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 131
            vbroadcastsd 4064(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12512(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4072(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4080(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4088(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 132
               vmovaps 12544(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 132
               vmovaps 12576(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 132
            vbroadcastsd 4096(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12608(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4104(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4112(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4120(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 133
               vmovaps 12640(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 133
               vmovaps 12672(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 133
            vbroadcastsd 4128(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12704(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4136(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4144(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4152(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 134
               vmovaps 12736(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 134
               vmovaps 12768(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 134
            vbroadcastsd 4160(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12800(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4168(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4176(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4184(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 135
               vmovaps 12832(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 135
               vmovaps 12864(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 135
            vbroadcastsd 4192(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12896(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4200(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4208(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4216(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 136
               vmovaps 12928(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 136
               vmovaps 12960(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 136
            vbroadcastsd 4224(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12992(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4232(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4240(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4248(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 137
               vmovaps 13024(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 137
               vmovaps 13056(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 137
            vbroadcastsd 4256(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13088(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4264(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4272(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4280(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 138
               vmovaps 13120(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 138
               vmovaps 13152(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 138
            vbroadcastsd 4288(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13184(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4296(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4304(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4312(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 139
               vmovaps 13216(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 139
               vmovaps 13248(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 139
            vbroadcastsd 4320(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13280(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4328(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4336(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4344(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 140
               vmovaps 13312(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 140
               vmovaps 13344(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 140
            vbroadcastsd 4352(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13376(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4360(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4368(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4376(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 141
               vmovaps 13408(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 141
               vmovaps 13440(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 141
            vbroadcastsd 4384(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13472(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4392(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4400(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4408(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 142
               vmovaps 13504(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 142
               vmovaps 13536(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 142
            vbroadcastsd 4416(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13568(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4424(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4432(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4440(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 143
               vmovaps 13600(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 143
               vmovaps 13632(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 143
            vbroadcastsd 4448(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13664(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4456(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4464(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4472(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 144
               vmovaps 13696(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 144
               vmovaps 13728(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 144
            vbroadcastsd 4480(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13760(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4488(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4496(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4504(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 145
               vmovaps 13792(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 145
               vmovaps 13824(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 145
            vbroadcastsd 4512(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13856(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4520(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4528(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4536(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 146
               vmovaps 13888(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 146
               vmovaps 13920(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 146
            vbroadcastsd 4544(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13952(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4552(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4560(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4568(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 147
               vmovaps 13984(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 147
               vmovaps 14016(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 147
            vbroadcastsd 4576(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14048(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4584(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4592(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4600(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 148
               vmovaps 14080(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 148
               vmovaps 14112(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 148
            vbroadcastsd 4608(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14144(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4616(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4624(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4632(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 149
               vmovaps 14176(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 149
               vmovaps 14208(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 149
            vbroadcastsd 4640(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14240(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4648(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4656(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4664(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 150
               vmovaps 14272(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 150
               vmovaps 14304(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 150
            vbroadcastsd 4672(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14336(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4680(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4688(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4696(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 151
               vmovaps 14368(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 151
               vmovaps 14400(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 151
            vbroadcastsd 4704(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14432(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4712(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4720(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4728(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 152
               vmovaps 14464(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 152
               vmovaps 14496(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 152
            vbroadcastsd 4736(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14528(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4744(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4752(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4760(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 153
               vmovaps 14560(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 153
               vmovaps 14592(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 153
            vbroadcastsd 4768(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14624(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4776(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4784(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4792(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 154
               vmovaps 14656(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 154
               vmovaps 14688(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 154
            vbroadcastsd 4800(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14720(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4808(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4816(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4824(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 155
               vmovaps 14752(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 155
               vmovaps 14784(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 155
            vbroadcastsd 4832(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14816(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4840(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4848(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4856(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 156
               vmovaps 14848(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 156
               vmovaps 14880(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 156
            vbroadcastsd 4864(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14912(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4872(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4880(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4888(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 157
               vmovaps 14944(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 157
               vmovaps 14976(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 157
            vbroadcastsd 4896(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15008(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4904(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4912(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4920(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 158
               vmovaps 15040(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 158
               vmovaps 15072(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 158
            vbroadcastsd 4928(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15104(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4936(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4944(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4952(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 159
               vmovaps 15136(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 159
               vmovaps 15168(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 159
            vbroadcastsd 4960(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15200(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 4968(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 4976(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 4984(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 160
               vmovaps 15232(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 160
               vmovaps 15264(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 160
            vbroadcastsd 4992(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15296(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5000(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5008(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5016(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 161
               vmovaps 15328(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 161
               vmovaps 15360(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 161
            vbroadcastsd 5024(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15392(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5032(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5040(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5048(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 162
               vmovaps 15424(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 162
               vmovaps 15456(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 162
            vbroadcastsd 5056(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15488(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5064(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5072(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5080(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 163
               vmovaps 15520(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 163
               vmovaps 15552(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 163
            vbroadcastsd 5088(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15584(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5096(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5104(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5112(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 164
               vmovaps 15616(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 164
               vmovaps 15648(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 164
            vbroadcastsd 5120(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15680(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5128(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5136(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5144(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 165
               vmovaps 15712(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 165
               vmovaps 15744(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 165
            vbroadcastsd 5152(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15776(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5160(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5168(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5176(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 166
               vmovaps 15808(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 166
               vmovaps 15840(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 166
            vbroadcastsd 5184(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15872(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5192(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5200(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5208(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 167
               vmovaps 15904(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 167
               vmovaps 15936(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 167
            vbroadcastsd 5216(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15968(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5224(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5232(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5240(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 168
               vmovaps 16000(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 168
               vmovaps 16032(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 168
            vbroadcastsd 5248(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16064(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5256(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5264(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5272(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 169
               vmovaps 16096(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 169
               vmovaps 16128(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 169
            vbroadcastsd 5280(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16160(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5288(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5296(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5304(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 170
               vmovaps 16192(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 170
               vmovaps 16224(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 170
            vbroadcastsd 5312(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16256(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5320(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5328(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5336(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 171
               vmovaps 16288(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 171
               vmovaps 16320(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 171
            vbroadcastsd 5344(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16352(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5352(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5360(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5368(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 172
               vmovaps 16384(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 172
               vmovaps 16416(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 172
            vbroadcastsd 5376(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16448(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5384(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5392(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5400(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 173
               vmovaps 16480(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 173
               vmovaps 16512(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 173
            vbroadcastsd 5408(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16544(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5416(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5424(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5432(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 174
               vmovaps 16576(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 174
               vmovaps 16608(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 174
            vbroadcastsd 5440(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16640(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5448(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5456(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5464(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 175
               vmovaps 16672(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 175
               vmovaps 16704(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 175
            vbroadcastsd 5472(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16736(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5480(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5488(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5496(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 176
               vmovaps 16768(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 176
               vmovaps 16800(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 176
            vbroadcastsd 5504(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16832(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5512(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5520(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5528(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 177
               vmovaps 16864(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 177
               vmovaps 16896(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 177
            vbroadcastsd 5536(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16928(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5544(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5552(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5560(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 178
               vmovaps 16960(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 178
               vmovaps 16992(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 178
            vbroadcastsd 5568(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17024(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5576(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5584(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5592(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 179
               vmovaps 17056(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 179
               vmovaps 17088(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 179
            vbroadcastsd 5600(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17120(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5608(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5616(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5624(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 180
               vmovaps 17152(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 180
               vmovaps 17184(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 180
            vbroadcastsd 5632(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17216(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5640(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5648(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5656(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 181
               vmovaps 17248(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 181
               vmovaps 17280(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 181
            vbroadcastsd 5664(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17312(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5672(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5680(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5688(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 182
               vmovaps 17344(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 182
               vmovaps 17376(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 182
            vbroadcastsd 5696(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17408(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5704(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5712(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5720(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 183
               vmovaps 17440(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 183
               vmovaps 17472(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 183
            vbroadcastsd 5728(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17504(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5736(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5744(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5752(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 184
               vmovaps 17536(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 184
               vmovaps 17568(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 184
            vbroadcastsd 5760(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17600(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5768(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5776(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5784(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 185
               vmovaps 17632(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 185
               vmovaps 17664(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 185
            vbroadcastsd 5792(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17696(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5800(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5808(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5816(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 186
               vmovaps 17728(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 186
               vmovaps 17760(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 186
            vbroadcastsd 5824(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17792(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5832(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5840(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5848(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 187
               vmovaps 17824(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 187
               vmovaps 17856(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 187
            vbroadcastsd 5856(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17888(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5864(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5872(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5880(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 188
               vmovaps 17920(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 188
               vmovaps 17952(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 188
            vbroadcastsd 5888(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17984(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5896(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5904(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5912(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 189
               vmovaps 18016(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 189
               vmovaps 18048(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 189
            vbroadcastsd 5920(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18080(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5928(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5936(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5944(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 190
               vmovaps 18112(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 190
               vmovaps 18144(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 190
            vbroadcastsd 5952(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18176(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5960(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 5968(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 5976(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 191
               vmovaps 18208(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 191
               vmovaps 18240(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 191
            vbroadcastsd 5984(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18272(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 5992(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6000(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6008(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 192
               vmovaps 18304(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 192
               vmovaps 18336(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 192
            vbroadcastsd 6016(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18368(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6024(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6032(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6040(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 193
               vmovaps 18400(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 193
               vmovaps 18432(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 193
            vbroadcastsd 6048(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18464(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6056(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6064(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6072(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 194
               vmovaps 18496(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 194
               vmovaps 18528(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 194
            vbroadcastsd 6080(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18560(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6088(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6096(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6104(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 195
               vmovaps 18592(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 195
               vmovaps 18624(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 195
            vbroadcastsd 6112(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18656(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6120(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6128(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6136(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 196
               vmovaps 18688(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 196
               vmovaps 18720(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 196
            vbroadcastsd 6144(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18752(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6152(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6160(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6168(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 197
               vmovaps 18784(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 197
               vmovaps 18816(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 197
            vbroadcastsd 6176(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18848(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6184(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6192(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6200(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 198
               vmovaps 18880(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 198
               vmovaps 18912(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 198
            vbroadcastsd 6208(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18944(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6216(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6224(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6232(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 199
               vmovaps 18976(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 199
               vmovaps 19008(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 199
            vbroadcastsd 6240(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19040(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6248(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6256(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6264(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 200
               vmovaps 19072(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 200
               vmovaps 19104(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 200
            vbroadcastsd 6272(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19136(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6280(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6288(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6296(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 201
               vmovaps 19168(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 201
               vmovaps 19200(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 201
            vbroadcastsd 6304(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19232(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6312(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6320(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6328(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 202
               vmovaps 19264(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 202
               vmovaps 19296(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 202
            vbroadcastsd 6336(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19328(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6344(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6352(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6360(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 203
               vmovaps 19360(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 203
               vmovaps 19392(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 203
            vbroadcastsd 6368(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19424(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6376(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6384(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6392(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 204
               vmovaps 19456(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 204
               vmovaps 19488(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 204
            vbroadcastsd 6400(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19520(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6408(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6416(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6424(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 205
               vmovaps 19552(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 205
               vmovaps 19584(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 205
            vbroadcastsd 6432(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19616(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6440(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6448(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6456(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 206
               vmovaps 19648(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 206
               vmovaps 19680(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 206
            vbroadcastsd 6464(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19712(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6472(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6480(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6488(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 207
               vmovaps 19744(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 207
               vmovaps 19776(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 207
            vbroadcastsd 6496(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19808(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6504(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6512(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6520(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 208
               vmovaps 19840(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 208
               vmovaps 19872(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 208
            vbroadcastsd 6528(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19904(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6536(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6544(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6552(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 209
               vmovaps 19936(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 209
               vmovaps 19968(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 209
            vbroadcastsd 6560(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20000(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6568(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6576(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6584(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 210
               vmovaps 20032(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 210
               vmovaps 20064(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 210
            vbroadcastsd 6592(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20096(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6600(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6608(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6616(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 211
               vmovaps 20128(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 211
               vmovaps 20160(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 211
            vbroadcastsd 6624(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20192(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6632(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6640(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6648(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 212
               vmovaps 20224(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 212
               vmovaps 20256(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 212
            vbroadcastsd 6656(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20288(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6664(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6672(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6680(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 213
               vmovaps 20320(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 213
               vmovaps 20352(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 213
            vbroadcastsd 6688(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20384(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6696(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6704(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6712(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 214
               vmovaps 20416(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 214
               vmovaps 20448(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 214
            vbroadcastsd 6720(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20480(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6728(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6736(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6744(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 215
               vmovaps 20512(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 215
               vmovaps 20544(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 215
            vbroadcastsd 6752(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20576(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6760(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6768(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6776(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 216
               vmovaps 20608(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 216
               vmovaps 20640(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 216
            vbroadcastsd 6784(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20672(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6792(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6800(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6808(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 217
               vmovaps 20704(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 217
               vmovaps 20736(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 217
            vbroadcastsd 6816(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20768(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6824(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6832(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6840(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 218
               vmovaps 20800(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 218
               vmovaps 20832(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 218
            vbroadcastsd 6848(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20864(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6856(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6864(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6872(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 219
               vmovaps 20896(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 219
               vmovaps 20928(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 219
            vbroadcastsd 6880(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20960(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6888(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6896(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6904(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 220
               vmovaps 20992(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 220
               vmovaps 21024(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 220
            vbroadcastsd 6912(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21056(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6920(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6928(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6936(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 221
               vmovaps 21088(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 221
               vmovaps 21120(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 221
            vbroadcastsd 6944(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21152(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6952(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6960(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 6968(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 222
               vmovaps 21184(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 222
               vmovaps 21216(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 222
            vbroadcastsd 6976(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21248(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 6984(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 6992(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7000(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 223
               vmovaps 21280(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 223
               vmovaps 21312(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 223
            vbroadcastsd 7008(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21344(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7016(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7024(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7032(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 224
               vmovaps 21376(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 224
               vmovaps 21408(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 224
            vbroadcastsd 7040(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21440(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7048(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7056(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7064(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 225
               vmovaps 21472(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 225
               vmovaps 21504(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 225
            vbroadcastsd 7072(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21536(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7080(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7088(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7096(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 226
               vmovaps 21568(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 226
               vmovaps 21600(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 226
            vbroadcastsd 7104(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21632(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7112(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7120(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7128(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 227
               vmovaps 21664(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 227
               vmovaps 21696(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 227
            vbroadcastsd 7136(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21728(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7144(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7152(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7160(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 228
               vmovaps 21760(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 228
               vmovaps 21792(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 228
            vbroadcastsd 7168(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21824(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7176(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7184(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7192(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 229
               vmovaps 21856(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 229
               vmovaps 21888(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 229
            vbroadcastsd 7200(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21920(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7208(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7216(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7224(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 230
               vmovaps 21952(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 230
               vmovaps 21984(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 230
            vbroadcastsd 7232(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22016(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7240(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7248(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7256(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 231
               vmovaps 22048(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 231
               vmovaps 22080(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 231
            vbroadcastsd 7264(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22112(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7272(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7280(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7288(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 232
               vmovaps 22144(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 232
               vmovaps 22176(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 232
            vbroadcastsd 7296(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22208(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7304(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7312(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7320(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 233
               vmovaps 22240(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 233
               vmovaps 22272(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 233
            vbroadcastsd 7328(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22304(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7336(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7344(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7352(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 234
               vmovaps 22336(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 234
               vmovaps 22368(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 234
            vbroadcastsd 7360(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22400(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7368(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7376(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7384(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 235
               vmovaps 22432(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 235
               vmovaps 22464(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 235
            vbroadcastsd 7392(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22496(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7400(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7408(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7416(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 236
               vmovaps 22528(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 236
               vmovaps 22560(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 236
            vbroadcastsd 7424(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22592(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7432(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7440(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7448(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 237
               vmovaps 22624(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 237
               vmovaps 22656(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 237
            vbroadcastsd 7456(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22688(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7464(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7472(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7480(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 238
               vmovaps 22720(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 238
               vmovaps 22752(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 238
            vbroadcastsd 7488(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22784(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7496(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7504(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7512(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 239
               vmovaps 22816(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 239
               vmovaps 22848(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 239
            vbroadcastsd 7520(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22880(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7528(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7536(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7544(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 240
               vmovaps 22912(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 240
               vmovaps 22944(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 240
            vbroadcastsd 7552(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22976(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7560(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7568(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7576(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 241
               vmovaps 23008(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 241
               vmovaps 23040(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 241
            vbroadcastsd 7584(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23072(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7592(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7600(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7608(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 242
               vmovaps 23104(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 242
               vmovaps 23136(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 242
            vbroadcastsd 7616(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23168(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7624(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7632(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7640(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 243
               vmovaps 23200(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 243
               vmovaps 23232(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 243
            vbroadcastsd 7648(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23264(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7656(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7664(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7672(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 244
               vmovaps 23296(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 244
               vmovaps 23328(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 244
            vbroadcastsd 7680(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23360(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7688(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7696(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7704(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 245
               vmovaps 23392(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 245
               vmovaps 23424(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 245
            vbroadcastsd 7712(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23456(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7720(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7728(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7736(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 246
               vmovaps 23488(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 246
               vmovaps 23520(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 246
            vbroadcastsd 7744(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23552(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7752(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7760(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7768(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 247
               vmovaps 23584(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 247
               vmovaps 23616(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 247
            vbroadcastsd 7776(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23648(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7784(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7792(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7800(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 248
               vmovaps 23680(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 248
               vmovaps 23712(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 248
            vbroadcastsd 7808(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23744(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7816(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7824(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7832(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 249
               vmovaps 23776(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 249
               vmovaps 23808(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 249
            vbroadcastsd 7840(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23840(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7848(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7856(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7864(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 250
               vmovaps 23872(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 250
               vmovaps 23904(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 250
            vbroadcastsd 7872(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23936(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7880(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7888(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7896(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 251
               vmovaps 23968(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 251
               vmovaps 24000(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 251
            vbroadcastsd 7904(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24032(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7912(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7920(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7928(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 252
               vmovaps 24064(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 252
               vmovaps 24096(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 252
            vbroadcastsd 7936(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24128(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7944(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7952(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7960(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 253
               vmovaps 24160(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 253
               vmovaps 24192(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 253
            vbroadcastsd 7968(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24224(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 7976(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 7984(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 7992(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 254
               vmovaps 24256(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 254
               vmovaps 24288(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 254
            vbroadcastsd 8000(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24320(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 8008(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 8016(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 8024(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 255
               vmovaps 24352(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 255
               vmovaps 24384(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 255
            vbroadcastsd 8032(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24416(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastsd 8040(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastsd 8048(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastsd 8056(pB), rB0
            FMAC rA0, rB0, rC03
            #if KB > 256
               vmovaps 24448(pA), rA0
            #endif
            FMAC rA1, rB0, rC13
            #if KB > 256
               vmovaps 24480(pA), rA1
            #endif
            FMAC rA2, rB0, rC23
         #endif

         add incPF, pfA
         #ifdef ATL_ARCH_AMDDRIVER
            #define USESTORE16 1
         #else
            #define USESTORE16 0
         #endif
         #ifndef BETA0
            VCOP -128(pC), rC00, rC00
         #endif
         #if USESTORE16
            movaps rC00x, -128(pC)
            vextractf128 $1, rC00, -112(pC)
         #else
            vmovapd rC00, -128(pC)
         #endif
         #ifndef BETA0
            VCOP -96(pC), rC10, rC10
         #endif
         #if USESTORE16
            movaps rC10x, -96(pC)
            vextractf128 $1, rC10, -80(pC)
         #else
            vmovapd rC10, -96(pC)
         #endif
         #ifndef BETA0
            VCOP -64(pC), rC20, rC20
         #endif
         #if USESTORE16
            movaps rC20x, -64(pC)
            vextractf128 $1, rC20, -48(pC)
         #else
            vmovapd rC20, -64(pC)
         #endif
         #ifndef BETA0
            VCOP -32(pC), rC01, rC01
         #endif
         #if USESTORE16
            movaps rC01x, -32(pC)
            vextractf128 $1, rC01, -16(pC)
         #else
            vmovapd rC01, -32(pC)
         #endif
         #ifndef BETA0
            VCOP (pC), rC11, rC11
         #endif
         #if USESTORE16
            movaps rC11x, (pC)
            vextractf128 $1, rC11, 16(pC)
         #else
            vmovapd rC11, (pC)
         #endif
         #ifndef BETA0
            VCOP 32(pC), rC21, rC21
         #endif
         #if USESTORE16
            movaps rC21x, 32(pC)
            vextractf128 $1, rC21, 48(pC)
         #else
            vmovapd rC21, 32(pC)
         #endif
         #ifndef BETA0
            VCOP 64(pC), rC02, rC02
         #endif
         #if USESTORE16
            movaps rC02x, 64(pC)
            vextractf128 $1, rC02, 80(pC)
         #else
            vmovapd rC02, 64(pC)
         #endif
         #ifndef BETA0
            VCOP 96(pC), rC12, rC12
         #endif
         #if USESTORE16
            movaps rC12x, 96(pC)
            vextractf128 $1, rC12, 112(pC)
         #else
            vmovapd rC12, 96(pC)
         #endif
         #ifndef BETA0
            VCOP 128(pC), rC22, rC22
         #endif
         #if USESTORE16
            movaps rC22x, 128(pC)
            vextractf128 $1, rC22, 144(pC)
         #else
            vmovapd rC22, 128(pC)
         #endif
         #ifndef BETA0
            VCOP 160(pC), rC03, rC03
         #endif
         #if USESTORE16
            movaps rC03x, 160(pC)
            vextractf128 $1, rC03, 176(pC)
         #else
            vmovapd rC03, 160(pC)
         #endif
         #ifndef BETA0
            VCOP 192(pC), rC13, rC13
         #endif
         #if USESTORE16
            movaps rC13x, 192(pC)
            vextractf128 $1, rC13, 208(pC)
         #else
            vmovapd rC13, 192(pC)
         #endif
         #ifndef BETA0
            VCOP 224(pC), rC23, rC23
         #endif
         #if USESTORE16
            movaps rC23x, 224(pC)
            vextractf128 $1, rC23, 240(pC)
         #else
            vmovapd rC23, 224(pC)
         #endif
         add $KB*4*8, pB     /* pB += KB*NU*sizeof */
         add $384, pC
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
