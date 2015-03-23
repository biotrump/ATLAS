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
#define FMAC vfmadd231ps   /* FMAC m256/r256, rs1, rd */
#if defined(BETAN) || defined(BETAn)
   #define BETAN1
#endif
#ifdef BETAN1
   #define VCOP vsubps
#else
   #define VCOP vaddps
#endif
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
   mov $192, r192
   mov pA, pA0
   mov pB, pB0
/*
 * incAm = 24*sizeof*K = 24*4*K = 3*4*8*K = 3*K*32
 */
   lea (KK, KK,2), incAm   /* incAm = 3*K */
   shl $5, incAm           /* incAm = 32*3*K */

   ALIGN16
   MLOOP:
      NLOOP:
/*
 *       First iteration peeled to handle init of rC
 */
         vbroadcastss      -128(pB), rB0
         vmovaps -128(pA), rA0
         vmulps rA0, rB0, rC00
         vmovaps -96(pA), rA1
         vmulps rA1, rB0, rC10
         vmovaps -64(pA), rA2
         vmulps rA2, rB0, rC20

         vbroadcastss      -124(pB), rB0
         vmulps rA0, rB0, rC01
            prefC(-128(pC))
         vmulps rA1, rB0, rC11
         vmulps rA2, rB0, rC21

         vbroadcastss      -120(pB), rB0
         vmulps rA0, rB0, rC02
         vmulps rA1, rB0, rC12
            prefC((pC))
         vmulps rA2, rB0, rC22

         vbroadcastss      -116(pB), rB0
         vmulps rA0, rB0, rC03
            vmovaps -32(pA), rA0
         vmulps rA1, rB0, rC13
            vmovaps (pA), rA1
         vmulps rA2, rB0, rC23

/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 1
            vbroadcastss -112(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 32(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss -108(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss -104(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss -100(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 64(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 96(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 2
            vbroadcastss -96(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 128(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss -92(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss -88(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss -84(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 160(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 192(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 3
            vbroadcastss -80(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 224(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss -76(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss -72(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss -68(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 256(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 288(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 4
            vbroadcastss -64(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 320(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss -60(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss -56(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss -52(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 352(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 384(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 5
            vbroadcastss -48(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 416(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss -44(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss -40(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss -36(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 448(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 480(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 6
            vbroadcastss -32(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 512(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss -28(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss -24(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss -20(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 544(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 576(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 7
            vbroadcastss -16(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 608(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss -12(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss -8(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss -4(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 640(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 672(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 8
            vbroadcastss 0(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 704(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 4(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 8(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 12(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 736(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 768(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 9
            vbroadcastss 16(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 800(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 20(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 24(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 28(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 832(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 864(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 10
            vbroadcastss 32(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 896(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 36(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 40(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 44(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 928(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 960(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 11
            vbroadcastss 48(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 992(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 52(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 56(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 60(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1024(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1056(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 12
            vbroadcastss 64(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1088(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 68(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 72(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 76(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1120(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1152(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 13
            vbroadcastss 80(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1184(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 84(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 88(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 92(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1216(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1248(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 14
            vbroadcastss 96(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1280(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 100(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 104(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 108(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1312(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1344(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 15
            vbroadcastss 112(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1376(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 116(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 120(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 124(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1408(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1440(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 16
            vbroadcastss 128(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1472(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 132(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 136(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 140(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1504(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1536(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 17
            vbroadcastss 144(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1568(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 148(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 152(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 156(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1600(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1632(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 18
            vbroadcastss 160(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1664(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 164(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 168(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 172(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1696(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1728(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 19
            vbroadcastss 176(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1760(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 180(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 184(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 188(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1792(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1824(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 20
            vbroadcastss 192(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1856(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 196(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 200(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 204(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1888(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 1920(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 21
            vbroadcastss 208(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 1952(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 212(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 216(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 220(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 1984(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2016(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 22
            vbroadcastss 224(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2048(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 228(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 232(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 236(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2080(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2112(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 23
            vbroadcastss 240(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2144(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 244(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 248(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 252(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2176(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2208(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 24
            vbroadcastss 256(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2240(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 260(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 264(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 268(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2272(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2304(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 25
            vbroadcastss 272(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2336(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 276(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 280(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 284(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2368(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2400(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 26
            vbroadcastss 288(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2432(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 292(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 296(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 300(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2464(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2496(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 27
            vbroadcastss 304(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2528(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 308(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 312(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 316(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2560(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2592(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 28
            vbroadcastss 320(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2624(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 324(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 328(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 332(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2656(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2688(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 29
            vbroadcastss 336(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2720(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 340(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 344(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 348(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2752(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2784(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 30
            vbroadcastss 352(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2816(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 356(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 360(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 364(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2848(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2880(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 31
            vbroadcastss 368(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 2912(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 372(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 376(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 380(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 2944(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 2976(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 32
            vbroadcastss 384(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3008(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 388(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 392(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 396(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3040(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3072(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 33
            vbroadcastss 400(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3104(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 404(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 408(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 412(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3136(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3168(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 34
            vbroadcastss 416(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3200(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 420(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 424(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 428(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3232(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3264(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 35
            vbroadcastss 432(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3296(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 436(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 440(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 444(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3328(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3360(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 36
            vbroadcastss 448(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3392(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 452(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 456(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 460(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3424(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3456(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 37
            vbroadcastss 464(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3488(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 468(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 472(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 476(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3520(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3552(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 38
            vbroadcastss 480(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3584(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 484(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 488(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 492(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3616(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3648(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 39
            vbroadcastss 496(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3680(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 500(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 504(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 508(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3712(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3744(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 40
            vbroadcastss 512(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3776(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 516(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 520(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 524(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3808(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3840(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 41
            vbroadcastss 528(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3872(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 532(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 536(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 540(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 3904(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 3936(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 42
            vbroadcastss 544(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 3968(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 548(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 552(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 556(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4000(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4032(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 43
            vbroadcastss 560(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4064(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 564(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 568(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 572(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4096(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4128(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 44
            vbroadcastss 576(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4160(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 580(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 584(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 588(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4192(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4224(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 45
            vbroadcastss 592(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4256(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 596(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 600(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 604(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4288(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4320(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 46
            vbroadcastss 608(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4352(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 612(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 616(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 620(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4384(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4416(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 47
            vbroadcastss 624(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4448(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 628(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 632(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 636(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4480(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4512(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 48
            vbroadcastss 640(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4544(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 644(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 648(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 652(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4576(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4608(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 49
            vbroadcastss 656(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4640(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 660(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 664(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 668(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4672(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4704(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 50
            vbroadcastss 672(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4736(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 676(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 680(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 684(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4768(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4800(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 51
            vbroadcastss 688(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4832(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 692(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 696(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 700(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4864(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4896(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 52
            vbroadcastss 704(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 4928(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 708(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 712(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 716(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 4960(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 4992(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 53
            vbroadcastss 720(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5024(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 724(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 728(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 732(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5056(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5088(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 54
            vbroadcastss 736(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5120(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 740(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 744(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 748(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5152(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5184(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 55
            vbroadcastss 752(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5216(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 756(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 760(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 764(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5248(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5280(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 56
            vbroadcastss 768(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5312(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 772(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 776(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 780(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5344(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5376(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 57
            vbroadcastss 784(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5408(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 788(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 792(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 796(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5440(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5472(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 58
            vbroadcastss 800(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5504(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 804(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 808(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 812(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5536(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5568(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 59
            vbroadcastss 816(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5600(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 820(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 824(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 828(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5632(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5664(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 60
            vbroadcastss 832(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5696(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 836(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 840(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 844(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5728(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5760(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 61
            vbroadcastss 848(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5792(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 852(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 856(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 860(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5824(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5856(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 62
            vbroadcastss 864(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5888(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 868(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 872(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 876(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 5920(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 5952(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 63
            vbroadcastss 880(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 5984(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 884(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 888(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 892(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6016(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6048(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 64
            vbroadcastss 896(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6080(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 900(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 904(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 908(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6112(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6144(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 65
            vbroadcastss 912(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6176(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 916(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 920(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 924(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6208(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6240(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 66
            vbroadcastss 928(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6272(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 932(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 936(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 940(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6304(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6336(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 67
            vbroadcastss 944(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6368(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 948(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 952(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 956(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6400(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6432(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 68
            vbroadcastss 960(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6464(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 964(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 968(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 972(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6496(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6528(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 69
            vbroadcastss 976(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6560(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 980(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 984(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 988(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6592(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6624(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 70
            vbroadcastss 992(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6656(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 996(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1000(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1004(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6688(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6720(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 71
            vbroadcastss 1008(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6752(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1012(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1016(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1020(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6784(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6816(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 72
            vbroadcastss 1024(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6848(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1028(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1032(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1036(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6880(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 6912(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 73
            vbroadcastss 1040(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 6944(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1044(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1048(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1052(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 6976(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7008(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 74
            vbroadcastss 1056(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7040(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1060(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1064(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1068(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7072(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7104(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 75
            vbroadcastss 1072(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7136(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1076(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1080(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1084(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7168(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7200(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 76
            vbroadcastss 1088(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7232(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1092(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1096(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1100(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7264(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7296(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 77
            vbroadcastss 1104(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7328(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1108(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1112(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1116(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7360(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7392(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 78
            vbroadcastss 1120(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7424(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1124(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1128(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1132(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7456(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7488(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 79
            vbroadcastss 1136(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7520(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1140(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1144(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1148(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7552(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7584(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 80
            vbroadcastss 1152(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7616(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1156(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1160(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1164(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7648(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7680(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 81
            vbroadcastss 1168(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7712(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1172(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1176(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1180(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7744(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7776(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 82
            vbroadcastss 1184(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7808(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1188(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1192(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1196(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7840(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7872(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 83
            vbroadcastss 1200(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 7904(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1204(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1208(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1212(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 7936(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 7968(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 84
            vbroadcastss 1216(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8000(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1220(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1224(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1228(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8032(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8064(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 85
            vbroadcastss 1232(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8096(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1236(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1240(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1244(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8128(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8160(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 86
            vbroadcastss 1248(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8192(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1252(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1256(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1260(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8224(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8256(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 87
            vbroadcastss 1264(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8288(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1268(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1272(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1276(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8320(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8352(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 88
            vbroadcastss 1280(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8384(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1284(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1288(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1292(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8416(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8448(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 89
            vbroadcastss 1296(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8480(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1300(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1304(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1308(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8512(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8544(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 90
            vbroadcastss 1312(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8576(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1316(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1320(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1324(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8608(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8640(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 91
            vbroadcastss 1328(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8672(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1332(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1336(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1340(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8704(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8736(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 92
            vbroadcastss 1344(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8768(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1348(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1352(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1356(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8800(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8832(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 93
            vbroadcastss 1360(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8864(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1364(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1368(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1372(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8896(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 8928(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 94
            vbroadcastss 1376(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 8960(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1380(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1384(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1388(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 8992(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9024(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 95
            vbroadcastss 1392(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9056(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1396(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1400(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1404(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9088(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9120(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 96
            vbroadcastss 1408(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9152(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1412(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1416(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1420(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9184(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9216(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 97
            vbroadcastss 1424(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9248(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1428(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1432(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1436(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9280(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9312(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 98
            vbroadcastss 1440(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9344(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1444(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1448(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1452(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9376(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9408(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 99
            vbroadcastss 1456(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9440(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1460(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1464(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1468(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9472(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9504(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 100
            vbroadcastss 1472(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9536(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1476(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1480(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1484(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9568(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9600(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 101
            vbroadcastss 1488(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9632(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1492(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1496(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1500(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9664(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9696(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 102
            vbroadcastss 1504(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9728(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1508(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1512(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1516(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9760(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9792(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 103
            vbroadcastss 1520(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9824(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1524(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1528(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1532(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9856(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9888(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 104
            vbroadcastss 1536(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 9920(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1540(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1544(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1548(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 9952(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 9984(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 105
            vbroadcastss 1552(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10016(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1556(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1560(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1564(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10048(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10080(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 106
            vbroadcastss 1568(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10112(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1572(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1576(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1580(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10144(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10176(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 107
            vbroadcastss 1584(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10208(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1588(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1592(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1596(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10240(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10272(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 108
            vbroadcastss 1600(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10304(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1604(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1608(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1612(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10336(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10368(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 109
            vbroadcastss 1616(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10400(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1620(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1624(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1628(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10432(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10464(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 110
            vbroadcastss 1632(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10496(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1636(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1640(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1644(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10528(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10560(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 111
            vbroadcastss 1648(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10592(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1652(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1656(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1660(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10624(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10656(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 112
            vbroadcastss 1664(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10688(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1668(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1672(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1676(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10720(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10752(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 113
            vbroadcastss 1680(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10784(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1684(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1688(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1692(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10816(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10848(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 114
            vbroadcastss 1696(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10880(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1700(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1704(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1708(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 10912(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 10944(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 115
            vbroadcastss 1712(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 10976(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1716(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1720(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1724(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11008(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11040(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 116
            vbroadcastss 1728(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11072(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1732(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1736(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1740(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11104(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11136(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 117
            vbroadcastss 1744(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11168(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1748(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1752(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1756(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11200(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11232(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 118
            vbroadcastss 1760(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11264(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1764(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1768(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1772(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11296(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11328(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 119
            vbroadcastss 1776(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11360(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1780(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1784(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1788(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11392(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11424(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 120
            vbroadcastss 1792(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11456(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1796(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1800(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1804(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11488(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11520(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 121
            vbroadcastss 1808(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11552(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1812(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1816(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1820(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11584(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11616(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 122
            vbroadcastss 1824(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11648(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1828(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1832(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1836(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11680(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11712(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 123
            vbroadcastss 1840(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11744(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1844(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1848(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1852(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11776(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11808(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 124
            vbroadcastss 1856(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11840(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1860(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1864(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1868(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11872(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 11904(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 125
            vbroadcastss 1872(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 11936(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1876(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1880(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1884(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 11968(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12000(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 126
            vbroadcastss 1888(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12032(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1892(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1896(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1900(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12064(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12096(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 127
            vbroadcastss 1904(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12128(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1908(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1912(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1916(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12160(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12192(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 128
            vbroadcastss 1920(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12224(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1924(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1928(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1932(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12256(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12288(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 129
            vbroadcastss 1936(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12320(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1940(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1944(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1948(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12352(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12384(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 130
            vbroadcastss 1952(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12416(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1956(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1960(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1964(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12448(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12480(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 131
            vbroadcastss 1968(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12512(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1972(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1976(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1980(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12544(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12576(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 132
            vbroadcastss 1984(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12608(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 1988(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 1992(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 1996(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12640(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12672(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 133
            vbroadcastss 2000(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12704(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2004(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2008(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2012(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12736(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12768(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 134
            vbroadcastss 2016(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12800(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2020(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2024(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2028(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12832(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12864(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 135
            vbroadcastss 2032(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12896(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2036(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2040(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2044(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 12928(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 12960(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 136
            vbroadcastss 2048(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 12992(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2052(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2056(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2060(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13024(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13056(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 137
            vbroadcastss 2064(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13088(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2068(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2072(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2076(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13120(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13152(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 138
            vbroadcastss 2080(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13184(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2084(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2088(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2092(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13216(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13248(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 139
            vbroadcastss 2096(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13280(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2100(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2104(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2108(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13312(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13344(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 140
            vbroadcastss 2112(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13376(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2116(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2120(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2124(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13408(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13440(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 141
            vbroadcastss 2128(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13472(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2132(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2136(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2140(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13504(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13536(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 142
            vbroadcastss 2144(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13568(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2148(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2152(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2156(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13600(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13632(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 143
            vbroadcastss 2160(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13664(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2164(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2168(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2172(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13696(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13728(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 144
            vbroadcastss 2176(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13760(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2180(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2184(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2188(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13792(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13824(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 145
            vbroadcastss 2192(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13856(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2196(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2200(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2204(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13888(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 13920(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 146
            vbroadcastss 2208(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 13952(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2212(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2216(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2220(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 13984(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14016(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 147
            vbroadcastss 2224(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14048(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2228(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2232(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2236(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14080(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14112(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 148
            vbroadcastss 2240(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14144(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2244(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2248(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2252(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14176(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14208(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 149
            vbroadcastss 2256(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14240(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2260(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2264(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2268(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14272(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14304(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 150
            vbroadcastss 2272(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14336(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2276(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2280(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2284(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14368(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14400(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 151
            vbroadcastss 2288(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14432(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2292(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2296(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2300(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14464(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14496(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 152
            vbroadcastss 2304(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14528(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2308(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2312(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2316(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14560(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14592(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 153
            vbroadcastss 2320(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14624(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2324(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2328(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2332(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14656(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14688(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 154
            vbroadcastss 2336(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14720(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2340(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2344(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2348(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14752(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14784(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 155
            vbroadcastss 2352(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14816(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2356(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2360(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2364(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14848(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14880(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 156
            vbroadcastss 2368(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 14912(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2372(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2376(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2380(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 14944(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 14976(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 157
            vbroadcastss 2384(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15008(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2388(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2392(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2396(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15040(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15072(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 158
            vbroadcastss 2400(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15104(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2404(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2408(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2412(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15136(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15168(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 159
            vbroadcastss 2416(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15200(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2420(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2424(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2428(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15232(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15264(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 160
            vbroadcastss 2432(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15296(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2436(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2440(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2444(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15328(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15360(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 161
            vbroadcastss 2448(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15392(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2452(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2456(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2460(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15424(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15456(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 162
            vbroadcastss 2464(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15488(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2468(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2472(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2476(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15520(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15552(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 163
            vbroadcastss 2480(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15584(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2484(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2488(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2492(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15616(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15648(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 164
            vbroadcastss 2496(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15680(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2500(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2504(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2508(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15712(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15744(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 165
            vbroadcastss 2512(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15776(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2516(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2520(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2524(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15808(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15840(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 166
            vbroadcastss 2528(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15872(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2532(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2536(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2540(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 15904(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 15936(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 167
            vbroadcastss 2544(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 15968(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2548(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2552(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2556(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16000(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16032(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 168
            vbroadcastss 2560(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16064(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2564(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2568(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2572(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16096(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16128(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 169
            vbroadcastss 2576(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16160(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2580(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2584(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2588(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16192(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16224(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 170
            vbroadcastss 2592(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16256(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2596(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2600(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2604(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16288(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16320(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 171
            vbroadcastss 2608(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16352(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2612(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2616(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2620(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16384(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16416(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 172
            vbroadcastss 2624(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16448(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2628(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2632(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2636(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16480(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16512(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 173
            vbroadcastss 2640(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16544(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2644(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2648(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2652(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16576(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16608(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 174
            vbroadcastss 2656(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16640(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2660(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2664(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2668(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16672(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16704(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 175
            vbroadcastss 2672(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16736(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2676(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2680(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2684(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16768(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16800(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 176
            vbroadcastss 2688(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16832(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2692(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2696(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2700(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16864(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16896(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 177
            vbroadcastss 2704(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 16928(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2708(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2712(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2716(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 16960(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 16992(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 178
            vbroadcastss 2720(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17024(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2724(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2728(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2732(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17056(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17088(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 179
            vbroadcastss 2736(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17120(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2740(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2744(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2748(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17152(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17184(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 180
            vbroadcastss 2752(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17216(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2756(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2760(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2764(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17248(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17280(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 181
            vbroadcastss 2768(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17312(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2772(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2776(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2780(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17344(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17376(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 182
            vbroadcastss 2784(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17408(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2788(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2792(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2796(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17440(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17472(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 183
            vbroadcastss 2800(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17504(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2804(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2808(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2812(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17536(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17568(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 184
            vbroadcastss 2816(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17600(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2820(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2824(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2828(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17632(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17664(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 185
            vbroadcastss 2832(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17696(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2836(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2840(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2844(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17728(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17760(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 186
            vbroadcastss 2848(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17792(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2852(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2856(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2860(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17824(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17856(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 187
            vbroadcastss 2864(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17888(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2868(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2872(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2876(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 17920(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 17952(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 188
            vbroadcastss 2880(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 17984(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2884(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2888(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2892(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18016(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18048(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 189
            vbroadcastss 2896(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18080(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2900(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2904(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2908(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18112(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18144(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 190
            vbroadcastss 2912(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18176(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2916(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2920(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2924(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18208(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18240(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 191
            vbroadcastss 2928(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18272(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2932(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2936(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2940(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18304(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18336(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 192
            vbroadcastss 2944(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18368(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2948(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2952(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2956(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18400(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18432(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 193
            vbroadcastss 2960(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18464(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2964(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2968(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2972(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18496(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18528(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 194
            vbroadcastss 2976(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18560(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2980(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 2984(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 2988(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18592(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18624(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 195
            vbroadcastss 2992(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18656(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 2996(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3000(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3004(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18688(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18720(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 196
            vbroadcastss 3008(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18752(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3012(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3016(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3020(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18784(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18816(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 197
            vbroadcastss 3024(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18848(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3028(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3032(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3036(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18880(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 18912(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 198
            vbroadcastss 3040(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 18944(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3044(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3048(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3052(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 18976(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19008(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 199
            vbroadcastss 3056(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19040(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3060(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3064(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3068(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19072(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19104(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 200
            vbroadcastss 3072(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19136(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3076(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3080(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3084(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19168(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19200(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 201
            vbroadcastss 3088(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19232(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3092(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3096(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3100(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19264(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19296(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 202
            vbroadcastss 3104(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19328(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3108(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3112(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3116(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19360(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19392(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 203
            vbroadcastss 3120(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19424(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3124(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3128(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3132(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19456(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19488(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 204
            vbroadcastss 3136(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19520(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3140(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3144(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3148(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19552(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19584(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 205
            vbroadcastss 3152(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19616(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3156(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3160(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3164(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19648(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19680(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 206
            vbroadcastss 3168(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19712(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3172(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3176(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3180(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19744(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19776(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 207
            vbroadcastss 3184(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19808(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3188(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3192(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3196(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19840(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19872(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 208
            vbroadcastss 3200(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 19904(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3204(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3208(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3212(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 19936(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 19968(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 209
            vbroadcastss 3216(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20000(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3220(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3224(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3228(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20032(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20064(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 210
            vbroadcastss 3232(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20096(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3236(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3240(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3244(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20128(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20160(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 211
            vbroadcastss 3248(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20192(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3252(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3256(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3260(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20224(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20256(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 212
            vbroadcastss 3264(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20288(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3268(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3272(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3276(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20320(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20352(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 213
            vbroadcastss 3280(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20384(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3284(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3288(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3292(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20416(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20448(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 214
            vbroadcastss 3296(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20480(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3300(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3304(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3308(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20512(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20544(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 215
            vbroadcastss 3312(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20576(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3316(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3320(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3324(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20608(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20640(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 216
            vbroadcastss 3328(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20672(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3332(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3336(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3340(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20704(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20736(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 217
            vbroadcastss 3344(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20768(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3348(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3352(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3356(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20800(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20832(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 218
            vbroadcastss 3360(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20864(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3364(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3368(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3372(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20896(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 20928(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 219
            vbroadcastss 3376(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 20960(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3380(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3384(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3388(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 20992(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21024(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 220
            vbroadcastss 3392(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21056(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3396(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3400(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3404(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21088(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21120(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 221
            vbroadcastss 3408(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21152(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3412(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3416(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3420(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21184(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21216(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 222
            vbroadcastss 3424(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21248(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3428(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3432(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3436(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21280(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21312(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 223
            vbroadcastss 3440(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21344(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3444(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3448(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3452(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21376(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21408(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 224
            vbroadcastss 3456(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21440(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3460(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3464(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3468(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21472(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21504(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 225
            vbroadcastss 3472(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21536(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3476(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3480(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3484(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21568(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21600(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 226
            vbroadcastss 3488(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21632(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3492(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3496(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3500(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21664(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21696(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 227
            vbroadcastss 3504(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21728(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3508(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3512(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3516(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21760(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21792(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 228
            vbroadcastss 3520(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21824(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3524(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3528(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3532(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21856(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21888(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 229
            vbroadcastss 3536(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 21920(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3540(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3544(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3548(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 21952(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 21984(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 230
            vbroadcastss 3552(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22016(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3556(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3560(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3564(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22048(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22080(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 231
            vbroadcastss 3568(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22112(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3572(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3576(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3580(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22144(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22176(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 232
            vbroadcastss 3584(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22208(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3588(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3592(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3596(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22240(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22272(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 233
            vbroadcastss 3600(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22304(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3604(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3608(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3612(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22336(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22368(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 234
            vbroadcastss 3616(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22400(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3620(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3624(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3628(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22432(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22464(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 235
            vbroadcastss 3632(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22496(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3636(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3640(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3644(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22528(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22560(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 236
            vbroadcastss 3648(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22592(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3652(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3656(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3660(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22624(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22656(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 237
            vbroadcastss 3664(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22688(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3668(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3672(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3676(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22720(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22752(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 238
            vbroadcastss 3680(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22784(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3684(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3688(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3692(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22816(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22848(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 239
            vbroadcastss 3696(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22880(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3700(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3704(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3708(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 22912(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 22944(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 240
            vbroadcastss 3712(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 22976(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3716(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3720(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3724(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23008(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23040(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 241
            vbroadcastss 3728(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23072(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3732(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3736(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3740(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23104(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23136(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 242
            vbroadcastss 3744(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23168(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3748(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3752(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3756(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23200(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23232(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 243
            vbroadcastss 3760(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23264(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3764(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3768(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3772(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23296(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23328(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 244
            vbroadcastss 3776(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23360(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3780(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3784(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3788(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23392(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23424(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 245
            vbroadcastss 3792(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23456(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3796(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3800(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3804(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23488(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23520(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 246
            vbroadcastss 3808(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23552(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3812(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3816(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3820(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23584(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23616(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 247
            vbroadcastss 3824(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23648(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3828(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3832(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3836(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23680(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23712(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 248
            vbroadcastss 3840(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23744(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3844(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3848(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3852(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23776(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23808(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 249
            vbroadcastss 3856(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23840(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3860(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3864(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3868(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23872(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 23904(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 250
            vbroadcastss 3872(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 23936(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3876(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3880(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3884(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 23968(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 24000(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 251
            vbroadcastss 3888(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24032(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3892(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3896(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3900(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 24064(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 24096(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 252
            vbroadcastss 3904(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24128(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3908(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3912(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3916(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 24160(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 24192(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 253
            vbroadcastss 3920(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24224(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3924(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3928(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3932(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 24256(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 24288(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 254
            vbroadcastss 3936(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24320(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3940(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3944(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3948(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 24352(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 24384(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 255
            vbroadcastss 3952(pB), rB0
            FMAC rA0, rB0, rC00
            vmovaps 24416(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            vbroadcastss 3956(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            vbroadcastss 3960(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            vbroadcastss 3964(pB), rB0
            FMAC rA0, rB0, rC03
               vmovaps 24448(pA), rA0
            FMAC rA1, rB0, rC13
               vmovaps 24480(pA), rA1
            FMAC rA2, rB0, rC23
         #endif

         #ifndef BETA0
            VCOP -128(pC), rC00, rC00
         #endif
         vmovaps rC00, -128(pC)
         #ifndef BETA0
            VCOP -96(pC), rC10, rC10
         #endif
         vmovaps rC10, -96(pC)
         #ifndef BETA0
            VCOP -64(pC), rC20, rC20
         #endif
         vmovaps rC20, -64(pC)
         #ifndef BETA0
            VCOP -32(pC), rC01, rC01
         #endif
         vmovaps rC01, -32(pC)
         #ifndef BETA0
            VCOP (pC), rC11, rC11
         #endif
         vmovaps rC11, (pC)
         #ifndef BETA0
            VCOP 32(pC), rC21, rC21
         #endif
         vmovaps rC21, 32(pC)
         #ifndef BETA0
            VCOP 64(pC), rC02, rC02
         #endif
         vmovaps rC02, 64(pC)
         #ifndef BETA0
            VCOP 96(pC), rC12, rC12
         #endif
         vmovaps rC12, 96(pC)
         #ifndef BETA0
            VCOP 128(pC), rC22, rC22
         #endif
         vmovaps rC22, 128(pC)
         #ifndef BETA0
            VCOP 160(pC), rC03, rC03
         #endif
         vmovaps rC03, 160(pC)
         #ifndef BETA0
            VCOP 192(pC), rC13, rC13
         #endif
         vmovaps rC13, 192(pC)
         #ifndef BETA0
            VCOP 224(pC), rC23, rC23
         #endif
         vmovaps rC23, 224(pC)
         add $KB*4*4, pB     /* pB += KB*NU*sizeof */
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
