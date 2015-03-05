#include "atlas_asm.h"

#define rB0     %xmm0
#define rA0     %xmm1
#define rA1     %xmm2
#define rA2     %xmm3
#define rC00    %xmm4
#define rC10    %xmm5
#define rC20    %xmm6
#define rC01    %xmm7
#define rC11    %xmm8
#define rC21    %xmm9
#define rC02    %xmm10
#define rC12    %xmm11
#define rC22    %xmm12
#define rC03    %xmm13
#define rC13    %xmm14
#define rC23    %xmm15
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
   #define VCOP subpd
#else
   #define VCOP addpd
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
   mov $192, r192
   mov pA, pA0
   mov pB, pB0
/*
 * incAm = 6*sizeof*K = 6*8*K = 3*2*8*K = 3*K*16
 */
   lea (KK, KK,2), incAm   /* incAm = 3*K */
   shl $4, incAm           /* incAm = 32*3*K */
   vxorps rA2, rA2, rA2

   ALIGN16
   MLOOP:
      NLOOP:
/*
 *       First iteration peeled to handle init of rC
 */
         vmovddup      -128(pB), rB0
         vmovaps -128(pA), rA0
         vmulpd rA0, rB0, rC00
         vmovaps -112(pA), rA1
         vmulpd rA1, rB0, rC10
         vmovaps -96(pA), rA2
         vmulpd rA2, rB0, rC20

         vmovddup      -120(pB), rB0
         vmulpd rA0, rB0, rC01
            prefC(-128(pC))
         vmulpd rA1, rB0, rC11
         vmulpd rA2, rB0, rC21

         vmovddup      -112(pB), rB0
         vmulpd rA0, rB0, rC02
         vmulpd rA1, rB0, rC12
            prefC((pC))
         vmulpd rA2, rB0, rC22

         vmovddup      -104(pB), rB0
         vmulpd rA0, rB0, rC03
            vmovapd -80(pA), rA0
         vmulpd rA1, rB0, rC13
            vmovapd -64(pA), rA1
         vmulpd rA2, rB0, rC23


/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 1
            movddup -96(pB), rB0
            FMAC rA0, rB0, rC00
            movapd -48(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup -88(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup -80(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup -72(pB), rB0
            FMAC rA0, rB0, rC03
               movapd -32(pA), rA0
            FMAC rA1, rB0, rC13
               movapd -16(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 2
            movddup -64(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 0(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup -56(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup -48(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup -40(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 16(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 32(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 3
            movddup -32(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 48(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup -24(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup -16(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup -8(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 64(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 80(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 4
            movddup 0(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 96(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 8(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 16(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 24(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 112(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 128(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 5
            movddup 32(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 144(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 40(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 48(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 56(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 160(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 176(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 6
            movddup 64(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 192(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 72(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 80(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 88(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 208(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 224(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 7
            movddup 96(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 240(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 104(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 112(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 120(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 256(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 272(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 8
            movddup 128(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 288(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 136(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 144(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 152(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 304(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 320(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 9
            movddup 160(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 336(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 168(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 176(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 184(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 352(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 368(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 10
            movddup 192(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 384(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 200(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 208(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 216(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 400(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 416(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 11
            movddup 224(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 432(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 232(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 240(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 248(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 448(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 464(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 12
            movddup 256(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 480(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 264(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 272(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 280(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 496(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 512(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 13
            movddup 288(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 528(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 296(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 304(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 312(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 544(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 560(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 14
            movddup 320(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 576(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 328(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 336(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 344(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 592(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 608(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 15
            movddup 352(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 624(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 360(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 368(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 376(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 640(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 656(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 16
            movddup 384(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 672(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 392(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 400(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 408(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 688(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 704(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 17
            movddup 416(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 720(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 424(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 432(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 440(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 736(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 752(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 18
            movddup 448(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 768(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 456(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 464(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 472(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 784(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 800(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 19
            movddup 480(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 816(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 488(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 496(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 504(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 832(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 848(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 20
            movddup 512(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 864(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 520(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 528(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 536(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 880(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 896(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 21
            movddup 544(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 912(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 552(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 560(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 568(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 928(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 944(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 22
            movddup 576(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 960(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 584(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 592(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 600(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 976(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 992(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 23
            movddup 608(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1008(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 616(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 624(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 632(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1024(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1040(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 24
            movddup 640(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1056(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 648(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 656(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 664(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1072(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1088(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 25
            movddup 672(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1104(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 680(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 688(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 696(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1120(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1136(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 26
            movddup 704(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1152(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 712(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 720(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 728(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1168(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1184(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 27
            movddup 736(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1200(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 744(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 752(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 760(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1216(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1232(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 28
            movddup 768(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1248(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 776(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 784(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 792(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1264(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1280(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 29
            movddup 800(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1296(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 808(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 816(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 824(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1312(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1328(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 30
            movddup 832(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1344(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 840(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 848(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 856(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1360(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1376(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 31
            movddup 864(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1392(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 872(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 880(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 888(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1408(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1424(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 32
            movddup 896(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1440(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 904(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 912(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 920(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1456(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1472(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 33
            movddup 928(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1488(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 936(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 944(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 952(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1504(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1520(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 34
            movddup 960(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1536(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 968(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 976(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 984(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1552(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1568(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 35
            movddup 992(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1584(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1000(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1008(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1016(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1600(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1616(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 36
            movddup 1024(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1632(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1032(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1040(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1048(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1648(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1664(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 37
            movddup 1056(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1680(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1064(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1072(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1080(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1696(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1712(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 38
            movddup 1088(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1728(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1096(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1104(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1112(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1744(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1760(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 39
            movddup 1120(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1776(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1128(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1136(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1144(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1792(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1808(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 40
            movddup 1152(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1824(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1160(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1168(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1176(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1840(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1856(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 41
            movddup 1184(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1872(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1192(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1200(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1208(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1888(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1904(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 42
            movddup 1216(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1920(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1224(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1232(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1240(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1936(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 1952(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 43
            movddup 1248(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 1968(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1256(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1264(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1272(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 1984(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2000(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 44
            movddup 1280(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2016(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1288(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1296(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1304(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2032(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2048(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 45
            movddup 1312(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2064(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1320(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1328(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1336(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2080(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2096(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 46
            movddup 1344(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2112(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1352(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1360(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1368(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2128(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2144(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 47
            movddup 1376(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2160(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1384(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1392(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1400(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2176(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2192(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 48
            movddup 1408(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2208(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1416(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1424(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1432(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2224(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2240(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 49
            movddup 1440(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2256(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1448(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1456(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1464(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2272(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2288(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 50
            movddup 1472(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2304(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1480(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1488(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1496(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2320(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2336(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 51
            movddup 1504(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2352(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1512(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1520(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1528(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2368(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2384(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 52
            movddup 1536(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2400(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1544(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1552(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1560(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2416(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2432(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 53
            movddup 1568(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2448(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1576(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1584(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1592(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2464(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2480(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 54
            movddup 1600(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2496(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1608(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1616(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1624(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2512(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2528(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 55
            movddup 1632(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2544(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1640(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1648(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1656(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2560(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2576(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 56
            movddup 1664(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2592(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1672(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1680(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1688(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2608(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2624(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 57
            movddup 1696(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2640(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1704(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1712(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1720(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2656(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2672(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 58
            movddup 1728(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2688(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1736(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1744(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1752(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2704(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2720(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 59
            movddup 1760(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2736(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1768(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1776(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1784(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2752(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2768(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 60
            movddup 1792(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2784(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1800(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1808(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1816(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2800(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2816(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 61
            movddup 1824(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2832(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1832(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1840(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1848(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2848(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2864(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 62
            movddup 1856(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2880(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1864(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1872(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1880(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2896(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2912(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 63
            movddup 1888(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2928(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1896(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1904(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1912(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2944(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 2960(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 64
            movddup 1920(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 2976(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1928(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1936(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1944(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 2992(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3008(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 65
            movddup 1952(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3024(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1960(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 1968(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 1976(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3040(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3056(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 66
            movddup 1984(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3072(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 1992(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2000(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2008(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3088(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3104(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 67
            movddup 2016(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3120(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2024(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2032(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2040(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3136(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3152(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 68
            movddup 2048(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3168(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2056(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2064(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2072(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3184(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3200(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 69
            movddup 2080(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3216(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2088(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2096(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2104(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3232(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3248(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 70
            movddup 2112(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3264(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2120(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2128(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2136(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3280(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3296(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 71
            movddup 2144(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3312(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2152(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2160(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2168(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3328(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3344(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 72
            movddup 2176(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3360(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2184(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2192(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2200(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3376(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3392(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 73
            movddup 2208(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3408(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2216(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2224(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2232(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3424(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3440(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 74
            movddup 2240(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3456(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2248(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2256(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2264(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3472(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3488(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 75
            movddup 2272(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3504(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2280(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2288(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2296(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3520(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3536(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 76
            movddup 2304(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3552(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2312(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2320(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2328(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3568(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3584(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 77
            movddup 2336(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3600(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2344(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2352(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2360(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3616(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3632(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 78
            movddup 2368(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3648(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2376(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2384(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2392(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3664(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3680(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 79
            movddup 2400(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3696(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2408(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2416(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2424(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3712(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3728(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 80
            movddup 2432(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3744(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2440(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2448(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2456(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3760(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3776(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 81
            movddup 2464(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3792(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2472(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2480(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2488(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3808(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3824(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 82
            movddup 2496(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3840(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2504(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2512(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2520(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3856(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3872(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 83
            movddup 2528(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3888(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2536(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2544(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2552(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3904(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3920(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 84
            movddup 2560(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3936(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2568(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2576(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2584(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 3952(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 3968(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 85
            movddup 2592(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 3984(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2600(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2608(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2616(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4000(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4016(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 86
            movddup 2624(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4032(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2632(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2640(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2648(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4048(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4064(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 87
            movddup 2656(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4080(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2664(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2672(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2680(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4096(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4112(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 88
            movddup 2688(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4128(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2696(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2704(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2712(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4144(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4160(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 89
            movddup 2720(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4176(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2728(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2736(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2744(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4192(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4208(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 90
            movddup 2752(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4224(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2760(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2768(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2776(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4240(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4256(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 91
            movddup 2784(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4272(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2792(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2800(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2808(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4288(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4304(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 92
            movddup 2816(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4320(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2824(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2832(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2840(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4336(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4352(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 93
            movddup 2848(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4368(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2856(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2864(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2872(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4384(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4400(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 94
            movddup 2880(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4416(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2888(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2896(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2904(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4432(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4448(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 95
            movddup 2912(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4464(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2920(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2928(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2936(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4480(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4496(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 96
            movddup 2944(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4512(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2952(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2960(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 2968(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4528(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4544(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 97
            movddup 2976(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4560(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 2984(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 2992(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3000(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4576(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4592(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 98
            movddup 3008(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4608(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3016(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3024(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3032(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4624(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4640(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 99
            movddup 3040(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4656(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3048(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3056(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3064(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4672(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4688(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 100
            movddup 3072(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4704(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3080(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3088(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3096(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4720(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4736(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 101
            movddup 3104(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4752(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3112(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3120(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3128(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4768(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4784(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 102
            movddup 3136(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4800(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3144(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3152(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3160(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4816(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4832(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 103
            movddup 3168(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4848(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3176(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3184(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3192(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4864(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4880(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 104
            movddup 3200(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4896(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3208(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3216(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3224(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4912(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4928(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 105
            movddup 3232(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4944(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3240(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3248(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3256(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 4960(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 4976(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 106
            movddup 3264(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 4992(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3272(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3280(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3288(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5008(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5024(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 107
            movddup 3296(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5040(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3304(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3312(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3320(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5056(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5072(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 108
            movddup 3328(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5088(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3336(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3344(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3352(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5104(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5120(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 109
            movddup 3360(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5136(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3368(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3376(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3384(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5152(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5168(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 110
            movddup 3392(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5184(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3400(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3408(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3416(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5200(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5216(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 111
            movddup 3424(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5232(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3432(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3440(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3448(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5248(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5264(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 112
            movddup 3456(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5280(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3464(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3472(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3480(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5296(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5312(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 113
            movddup 3488(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5328(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3496(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3504(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3512(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5344(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5360(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 114
            movddup 3520(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5376(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3528(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3536(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3544(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5392(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5408(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 115
            movddup 3552(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5424(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3560(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3568(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3576(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5440(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5456(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 116
            movddup 3584(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5472(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3592(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3600(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3608(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5488(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5504(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 117
            movddup 3616(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5520(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3624(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3632(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3640(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5536(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5552(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 118
            movddup 3648(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5568(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3656(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3664(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3672(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5584(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5600(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 119
            movddup 3680(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5616(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3688(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3696(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3704(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5632(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5648(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 120
            movddup 3712(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5664(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3720(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3728(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3736(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5680(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5696(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 121
            movddup 3744(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5712(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3752(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3760(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3768(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5728(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5744(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 122
            movddup 3776(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5760(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3784(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3792(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3800(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5776(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5792(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 123
            movddup 3808(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5808(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3816(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3824(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3832(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5824(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5840(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 124
            movddup 3840(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5856(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3848(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3856(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3864(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5872(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5888(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 125
            movddup 3872(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5904(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3880(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3888(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3896(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5920(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5936(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 126
            movddup 3904(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 5952(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3912(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3920(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3928(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 5968(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 5984(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 127
            movddup 3936(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6000(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3944(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3952(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3960(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6016(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6032(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 128
            movddup 3968(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6048(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 3976(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 3984(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 3992(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6064(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6080(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 129
            movddup 4000(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6096(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4008(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4016(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4024(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6112(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6128(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 130
            movddup 4032(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6144(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4040(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4048(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4056(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6160(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6176(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 131
            movddup 4064(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6192(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4072(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4080(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4088(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6208(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6224(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 132
            movddup 4096(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6240(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4104(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4112(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4120(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6256(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6272(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 133
            movddup 4128(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6288(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4136(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4144(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4152(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6304(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6320(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 134
            movddup 4160(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6336(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4168(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4176(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4184(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6352(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6368(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 135
            movddup 4192(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6384(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4200(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4208(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4216(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6400(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6416(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 136
            movddup 4224(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6432(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4232(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4240(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4248(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6448(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6464(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 137
            movddup 4256(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6480(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4264(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4272(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4280(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6496(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6512(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 138
            movddup 4288(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6528(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4296(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4304(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4312(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6544(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6560(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 139
            movddup 4320(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6576(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4328(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4336(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4344(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6592(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6608(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 140
            movddup 4352(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6624(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4360(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4368(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4376(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6640(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6656(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 141
            movddup 4384(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6672(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4392(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4400(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4408(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6688(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6704(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 142
            movddup 4416(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6720(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4424(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4432(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4440(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6736(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6752(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 143
            movddup 4448(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6768(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4456(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4464(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4472(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6784(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6800(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 144
            movddup 4480(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6816(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4488(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4496(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4504(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6832(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6848(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 145
            movddup 4512(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6864(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4520(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4528(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4536(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6880(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6896(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 146
            movddup 4544(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6912(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4552(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4560(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4568(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6928(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6944(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 147
            movddup 4576(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 6960(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4584(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4592(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4600(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 6976(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 6992(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 148
            movddup 4608(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7008(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4616(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4624(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4632(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7024(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7040(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 149
            movddup 4640(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7056(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4648(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4656(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4664(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7072(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7088(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 150
            movddup 4672(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7104(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4680(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4688(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4696(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7120(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7136(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 151
            movddup 4704(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7152(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4712(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4720(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4728(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7168(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7184(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 152
            movddup 4736(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7200(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4744(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4752(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4760(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7216(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7232(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 153
            movddup 4768(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7248(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4776(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4784(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4792(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7264(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7280(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 154
            movddup 4800(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7296(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4808(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4816(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4824(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7312(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7328(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 155
            movddup 4832(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7344(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4840(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4848(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4856(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7360(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7376(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 156
            movddup 4864(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7392(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4872(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4880(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4888(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7408(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7424(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 157
            movddup 4896(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7440(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4904(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4912(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4920(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7456(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7472(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 158
            movddup 4928(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7488(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4936(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4944(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4952(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7504(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7520(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 159
            movddup 4960(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7536(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 4968(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 4976(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 4984(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7552(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7568(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 160
            movddup 4992(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7584(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5000(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5008(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5016(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7600(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7616(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 161
            movddup 5024(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7632(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5032(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5040(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5048(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7648(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7664(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 162
            movddup 5056(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7680(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5064(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5072(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5080(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7696(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7712(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 163
            movddup 5088(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7728(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5096(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5104(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5112(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7744(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7760(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 164
            movddup 5120(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7776(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5128(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5136(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5144(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7792(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7808(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 165
            movddup 5152(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7824(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5160(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5168(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5176(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7840(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7856(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 166
            movddup 5184(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7872(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5192(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5200(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5208(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7888(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7904(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 167
            movddup 5216(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7920(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5224(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5232(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5240(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7936(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 7952(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 168
            movddup 5248(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 7968(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5256(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5264(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5272(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 7984(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8000(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 169
            movddup 5280(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8016(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5288(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5296(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5304(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8032(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8048(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 170
            movddup 5312(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8064(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5320(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5328(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5336(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8080(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8096(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 171
            movddup 5344(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8112(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5352(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5360(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5368(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8128(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8144(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 172
            movddup 5376(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8160(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5384(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5392(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5400(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8176(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8192(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 173
            movddup 5408(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8208(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5416(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5424(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5432(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8224(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8240(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 174
            movddup 5440(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8256(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5448(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5456(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5464(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8272(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8288(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 175
            movddup 5472(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8304(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5480(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5488(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5496(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8320(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8336(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 176
            movddup 5504(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8352(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5512(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5520(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5528(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8368(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8384(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 177
            movddup 5536(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8400(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5544(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5552(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5560(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8416(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8432(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 178
            movddup 5568(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8448(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5576(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5584(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5592(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8464(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8480(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 179
            movddup 5600(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8496(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5608(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5616(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5624(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8512(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8528(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 180
            movddup 5632(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8544(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5640(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5648(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5656(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8560(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8576(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 181
            movddup 5664(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8592(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5672(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5680(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5688(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8608(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8624(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 182
            movddup 5696(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8640(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5704(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5712(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5720(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8656(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8672(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 183
            movddup 5728(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8688(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5736(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5744(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5752(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8704(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8720(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 184
            movddup 5760(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8736(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5768(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5776(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5784(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8752(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8768(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 185
            movddup 5792(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8784(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5800(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5808(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5816(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8800(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8816(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 186
            movddup 5824(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8832(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5832(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5840(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5848(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8848(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8864(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 187
            movddup 5856(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8880(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5864(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5872(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5880(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8896(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8912(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 188
            movddup 5888(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8928(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5896(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5904(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5912(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8944(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 8960(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 189
            movddup 5920(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 8976(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5928(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5936(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5944(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 8992(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9008(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 190
            movddup 5952(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9024(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5960(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 5968(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 5976(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9040(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9056(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 191
            movddup 5984(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9072(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 5992(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6000(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6008(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9088(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9104(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 192
            movddup 6016(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9120(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6024(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6032(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6040(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9136(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9152(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 193
            movddup 6048(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9168(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6056(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6064(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6072(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9184(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9200(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 194
            movddup 6080(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9216(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6088(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6096(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6104(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9232(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9248(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 195
            movddup 6112(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9264(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6120(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6128(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6136(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9280(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9296(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 196
            movddup 6144(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9312(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6152(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6160(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6168(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9328(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9344(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 197
            movddup 6176(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9360(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6184(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6192(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6200(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9376(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9392(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 198
            movddup 6208(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9408(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6216(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6224(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6232(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9424(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9440(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 199
            movddup 6240(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9456(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6248(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6256(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6264(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9472(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9488(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 200
            movddup 6272(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9504(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6280(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6288(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6296(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9520(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9536(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 201
            movddup 6304(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9552(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6312(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6320(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6328(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9568(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9584(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 202
            movddup 6336(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9600(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6344(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6352(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6360(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9616(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9632(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 203
            movddup 6368(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9648(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6376(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6384(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6392(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9664(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9680(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 204
            movddup 6400(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9696(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6408(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6416(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6424(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9712(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9728(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 205
            movddup 6432(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9744(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6440(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6448(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6456(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9760(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9776(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 206
            movddup 6464(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9792(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6472(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6480(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6488(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9808(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9824(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 207
            movddup 6496(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9840(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6504(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6512(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6520(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9856(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9872(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 208
            movddup 6528(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9888(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6536(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6544(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6552(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9904(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9920(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 209
            movddup 6560(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9936(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6568(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6576(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6584(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 9952(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 9968(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 210
            movddup 6592(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 9984(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6600(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6608(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6616(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10000(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10016(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 211
            movddup 6624(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10032(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6632(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6640(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6648(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10048(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10064(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 212
            movddup 6656(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10080(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6664(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6672(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6680(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10096(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10112(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 213
            movddup 6688(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10128(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6696(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6704(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6712(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10144(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10160(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 214
            movddup 6720(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10176(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6728(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6736(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6744(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10192(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10208(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 215
            movddup 6752(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10224(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6760(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6768(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6776(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10240(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10256(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 216
            movddup 6784(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10272(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6792(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6800(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6808(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10288(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10304(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 217
            movddup 6816(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10320(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6824(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6832(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6840(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10336(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10352(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 218
            movddup 6848(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10368(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6856(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6864(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6872(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10384(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10400(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 219
            movddup 6880(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10416(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6888(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6896(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6904(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10432(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10448(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 220
            movddup 6912(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10464(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6920(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6928(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6936(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10480(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10496(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 221
            movddup 6944(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10512(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6952(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6960(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 6968(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10528(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10544(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 222
            movddup 6976(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10560(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 6984(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 6992(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7000(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10576(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10592(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 223
            movddup 7008(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10608(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7016(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7024(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7032(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10624(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10640(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 224
            movddup 7040(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10656(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7048(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7056(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7064(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10672(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10688(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 225
            movddup 7072(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10704(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7080(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7088(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7096(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10720(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10736(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 226
            movddup 7104(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10752(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7112(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7120(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7128(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10768(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10784(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 227
            movddup 7136(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10800(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7144(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7152(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7160(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10816(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10832(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 228
            movddup 7168(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10848(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7176(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7184(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7192(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10864(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10880(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 229
            movddup 7200(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10896(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7208(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7216(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7224(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10912(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10928(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 230
            movddup 7232(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10944(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7240(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7248(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7256(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 10960(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 10976(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 231
            movddup 7264(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 10992(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7272(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7280(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7288(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11008(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11024(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 232
            movddup 7296(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11040(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7304(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7312(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7320(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11056(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11072(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 233
            movddup 7328(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11088(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7336(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7344(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7352(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11104(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11120(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 234
            movddup 7360(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11136(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7368(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7376(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7384(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11152(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11168(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 235
            movddup 7392(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11184(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7400(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7408(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7416(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11200(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11216(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 236
            movddup 7424(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11232(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7432(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7440(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7448(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11248(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11264(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 237
            movddup 7456(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11280(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7464(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7472(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7480(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11296(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11312(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 238
            movddup 7488(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11328(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7496(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7504(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7512(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11344(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11360(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 239
            movddup 7520(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11376(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7528(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7536(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7544(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11392(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11408(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 240
            movddup 7552(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11424(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7560(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7568(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7576(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11440(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11456(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 241
            movddup 7584(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11472(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7592(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7600(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7608(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11488(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11504(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 242
            movddup 7616(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11520(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7624(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7632(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7640(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11536(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11552(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 243
            movddup 7648(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11568(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7656(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7664(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7672(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11584(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11600(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 244
            movddup 7680(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11616(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7688(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7696(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7704(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11632(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11648(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 245
            movddup 7712(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11664(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7720(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7728(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7736(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11680(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11696(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 246
            movddup 7744(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11712(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7752(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7760(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7768(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11728(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11744(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 247
            movddup 7776(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11760(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7784(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7792(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7800(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11776(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11792(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 248
            movddup 7808(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11808(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7816(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7824(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7832(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11824(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11840(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 249
            movddup 7840(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11856(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7848(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7856(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7864(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11872(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11888(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 250
            movddup 7872(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11904(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7880(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7888(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7896(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11920(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11936(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 251
            movddup 7904(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 11952(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7912(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7920(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7928(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 11968(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 11984(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 252
            movddup 7936(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 12000(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7944(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7952(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7960(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 12016(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 12032(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 253
            movddup 7968(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 12048(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 7976(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 7984(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 7992(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 12064(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 12080(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 254
            movddup 8000(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 12096(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 8008(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 8016(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 8024(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 12112(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 12128(pA), rA1
            FMAC rA2, rB0, rC23
         #endif
         #if KB > 255
            movddup 8032(pB), rB0
            FMAC rA0, rB0, rC00
            movapd 12144(pA), rA2
            FMAC rA1, rB0, rC10
            FMAC rA2, rB0, rC20

            movddup 8040(pB), rB0
            FMAC rA0, rB0, rC01
            FMAC rA1, rB0, rC11
            FMAC rA2, rB0, rC21

            movddup 8048(pB), rB0
            FMAC rA0, rB0, rC02
            FMAC rA1, rB0, rC12
            FMAC rA2, rB0, rC22

            movddup 8056(pB), rB0
            FMAC rA0, rB0, rC03
               movapd 12160(pA), rA0
            FMAC rA1, rB0, rC13
               movapd 12176(pA), rA1
            FMAC rA2, rB0, rC23
         #endif

         #ifndef BETA0
            VCOP -128(pC), rC00
         #endif
         movapd rC00, -128(pC)
         #ifndef BETA0
            VCOP -112(pC), rC10
         #endif
         movapd rC10, -112(pC)
         #ifndef BETA0
            VCOP -96(pC), rC20
         #endif
         movapd rC20, -96(pC)
         #ifndef BETA0
            VCOP -80(pC), rC01
         #endif
         movapd rC01, -80(pC)
         #ifndef BETA0
            VCOP -64(pC), rC11
         #endif
         movapd rC11, -64(pC)
         #ifndef BETA0
            VCOP -48(pC), rC21
         #endif
         movapd rC21, -48(pC)
         #ifndef BETA0
            VCOP -32(pC), rC02
         #endif
         movapd rC02, -32(pC)
         #ifndef BETA0
            VCOP -16(pC), rC12
         #endif
         movapd rC12, -16(pC)
         #ifndef BETA0
            VCOP (pC), rC22
         #endif
         movapd rC22, (pC)
         #ifndef BETA0
            VCOP 16(pC), rC03
         #endif
         movapd rC03, 16(pC)
         #ifndef BETA0
            VCOP 32(pC), rC13
         #endif
         movapd rC13, 32(pC)
         #ifndef BETA0
            VCOP 48(pC), rC23
         #endif
         movapd rC23, 48(pC)
         add $KB*4*8, pB     /* pB += KB*NU*sizeof */
         add $192, pC
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
