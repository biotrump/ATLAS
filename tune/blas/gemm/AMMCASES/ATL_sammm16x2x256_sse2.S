#include "atlas_asm.h"

#define rm0     %xmm0
#define rB0     %xmm1
#define rB1     %xmm2
#define rA0     %xmm3
#define rA1     %xmm4
#define rA2     %xmm5
#define rA3     %xmm6
#define rC00    %xmm7
#define rC10    %xmm8
#define rC20    %xmm9
#define rC30    %xmm10
#define rC01    %xmm11
#define rC11    %xmm12
#define rC21    %xmm13
#define rC31    %xmm14
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
   #define VCOP subps
#else
   #define VCOP addps
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
   cmp pfB, pB
   CMOVE pfA, pfB
   mov KK, KK0
   sub $-128, pfB
   sub $-128, pfA
   sub $-128, pA
   sub $-128, pB
   mov $24, r24
   mov pA, pA0
   mov pB, pB0
   mov $8*4*8, incPF
/*
 * incAm = MU*K*sizeof = 16*K*4 = 64*K
 */
   mov KK, incAm
   shl $6, incAm           /* incAm = 64*K */

   ALIGN16
   MLOOP:
      NLOOP:
/*
 *       First iteration peeled to handle init of rC
 */
         xorps rC00, rC00
         movaps rC00, rC10
         xorps rC20, rC20
         movaps rC00, rC30
         xorps rC01, rC01
         movaps rC00, rC11
         xorps rC21, rC21
         movaps rC00, rC31

         movddup -128(pB), rB0   /* b1 b0 b1 b0 */
         movaps -128(pA), rA0
         movaps -112(pA), rA1
         movaps -96(pA), rA2
         movaps -80(pA), rA3
/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         ALIGN16
         #if KB > 0
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 1
               movddup -120(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 1
               movaps -64(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 1
               movaps -48(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 1
               movaps -32(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 1
               movaps -16(pA), rA3
            #endif
         #endif
         #if KB > 1
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 2
               movddup -112(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 2
               movaps 0(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 2
               movaps 16(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 2
               movaps 32(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 2
               movaps 48(pA), rA3
            #endif
         #endif
         #if KB > 2
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 3
               movddup -104(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 3
               movaps 64(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 3
               movaps 80(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 3
               movaps 96(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 3
               movaps 112(pA), rA3
            #endif
         #endif
         #if KB > 3
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 4
               movddup -96(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 4
               movaps 128(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 4
               movaps 144(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 4
               movaps 160(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 4
               movaps 176(pA), rA3
            #endif
         #endif
         #if KB > 4
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 5
               movddup -88(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 5
               movaps 192(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 5
               movaps 208(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 5
               movaps 224(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 5
               movaps 240(pA), rA3
            #endif
         #endif
         #if KB > 5
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 6
               movddup -80(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 6
               movaps 256(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 6
               movaps 272(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 6
               movaps 288(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 6
               movaps 304(pA), rA3
            #endif
         #endif
         #if KB > 6
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 7
               movddup -72(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 7
               movaps 320(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 7
               movaps 336(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 7
               movaps 352(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 7
               movaps 368(pA), rA3
            #endif
         #endif
         #if KB > 7
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 8
               movddup -64(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 8
               movaps 384(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 8
               movaps 400(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 8
               movaps 416(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 8
               movaps 432(pA), rA3
            #endif
         #endif
         #if KB > 8
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 9
               movddup -56(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 9
               movaps 448(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 9
               movaps 464(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 9
               movaps 480(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 9
               movaps 496(pA), rA3
            #endif
         #endif
         #if KB > 9
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 10
               movddup -48(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 10
               movaps 512(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 10
               movaps 528(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 10
               movaps 544(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 10
               movaps 560(pA), rA3
            #endif
         #endif
         #if KB > 10
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 11
               movddup -40(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 11
               movaps 576(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 11
               movaps 592(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 11
               movaps 608(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 11
               movaps 624(pA), rA3
            #endif
         #endif
         #if KB > 11
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 12
               movddup -32(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 12
               movaps 640(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 12
               movaps 656(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 12
               movaps 672(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 12
               movaps 688(pA), rA3
            #endif
         #endif
         #if KB > 12
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 13
               movddup -24(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 13
               movaps 704(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 13
               movaps 720(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 13
               movaps 736(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 13
               movaps 752(pA), rA3
            #endif
         #endif
         #if KB > 13
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 14
               movddup -16(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 14
               movaps 768(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 14
               movaps 784(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 14
               movaps 800(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 14
               movaps 816(pA), rA3
            #endif
         #endif
         #if KB > 14
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 15
               movddup -8(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 15
               movaps 832(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 15
               movaps 848(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 15
               movaps 864(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 15
               movaps 880(pA), rA3
            #endif
         #endif
         #if KB > 15
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 16
               movddup 0(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 16
               movaps 896(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 16
               movaps 912(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 16
               movaps 928(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 16
               movaps 944(pA), rA3
            #endif
         #endif
         #if KB > 16
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 17
               movddup 8(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 17
               movaps 960(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 17
               movaps 976(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 17
               movaps 992(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 17
               movaps 1008(pA), rA3
            #endif
         #endif
         #if KB > 17
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 18
               movddup 16(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 18
               movaps 1024(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 18
               movaps 1040(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 18
               movaps 1056(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 18
               movaps 1072(pA), rA3
            #endif
         #endif
         #if KB > 18
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 19
               movddup 24(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 19
               movaps 1088(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 19
               movaps 1104(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 19
               movaps 1120(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 19
               movaps 1136(pA), rA3
            #endif
         #endif
         #if KB > 19
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 20
               movddup 32(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 20
               movaps 1152(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 20
               movaps 1168(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 20
               movaps 1184(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 20
               movaps 1200(pA), rA3
            #endif
         #endif
         #if KB > 20
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 21
               movddup 40(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 21
               movaps 1216(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 21
               movaps 1232(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 21
               movaps 1248(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 21
               movaps 1264(pA), rA3
            #endif
         #endif
         #if KB > 21
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 22
               movddup 48(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 22
               movaps 1280(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 22
               movaps 1296(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 22
               movaps 1312(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 22
               movaps 1328(pA), rA3
            #endif
         #endif
         #if KB > 22
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 23
               movddup 56(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 23
               movaps 1344(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 23
               movaps 1360(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 23
               movaps 1376(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 23
               movaps 1392(pA), rA3
            #endif
         #endif
         #if KB > 23
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 24
               movddup 64(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 24
               movaps 1408(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 24
               movaps 1424(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 24
               movaps 1440(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 24
               movaps 1456(pA), rA3
            #endif
         #endif
         #if KB > 24
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 25
               movddup 72(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 25
               movaps 1472(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 25
               movaps 1488(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 25
               movaps 1504(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 25
               movaps 1520(pA), rA3
            #endif
         #endif
         #if KB > 25
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 26
               movddup 80(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 26
               movaps 1536(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 26
               movaps 1552(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 26
               movaps 1568(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 26
               movaps 1584(pA), rA3
            #endif
         #endif
         #if KB > 26
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 27
               movddup 88(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 27
               movaps 1600(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 27
               movaps 1616(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 27
               movaps 1632(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 27
               movaps 1648(pA), rA3
            #endif
         #endif
         #if KB > 27
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 28
               movddup 96(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 28
               movaps 1664(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 28
               movaps 1680(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 28
               movaps 1696(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 28
               movaps 1712(pA), rA3
            #endif
         #endif
         #if KB > 28
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 29
               movddup 104(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 29
               movaps 1728(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 29
               movaps 1744(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 29
               movaps 1760(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 29
               movaps 1776(pA), rA3
            #endif
         #endif
         #if KB > 29
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 30
               movddup 112(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 30
               movaps 1792(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 30
               movaps 1808(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 30
               movaps 1824(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 30
               movaps 1840(pA), rA3
            #endif
         #endif
         #if KB > 30
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 31
               movddup 120(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 31
               movaps 1856(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 31
               movaps 1872(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 31
               movaps 1888(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 31
               movaps 1904(pA), rA3
            #endif
         #endif
         #if KB > 31
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 32
               movddup 128(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 32
               movaps 1920(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 32
               movaps 1936(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 32
               movaps 1952(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 32
               movaps 1968(pA), rA3
            #endif
         #endif
         #if KB > 32
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 33
               movddup 136(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 33
               movaps 1984(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 33
               movaps 2000(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 33
               movaps 2016(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 33
               movaps 2032(pA), rA3
            #endif
         #endif
         #if KB > 33
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 34
               movddup 144(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 34
               movaps 2048(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 34
               movaps 2064(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 34
               movaps 2080(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 34
               movaps 2096(pA), rA3
            #endif
         #endif
         #if KB > 34
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 35
               movddup 152(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 35
               movaps 2112(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 35
               movaps 2128(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 35
               movaps 2144(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 35
               movaps 2160(pA), rA3
            #endif
         #endif
         #if KB > 35
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 36
               movddup 160(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 36
               movaps 2176(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 36
               movaps 2192(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 36
               movaps 2208(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 36
               movaps 2224(pA), rA3
            #endif
         #endif
         #if KB > 36
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 37
               movddup 168(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 37
               movaps 2240(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 37
               movaps 2256(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 37
               movaps 2272(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 37
               movaps 2288(pA), rA3
            #endif
         #endif
         #if KB > 37
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 38
               movddup 176(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 38
               movaps 2304(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 38
               movaps 2320(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 38
               movaps 2336(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 38
               movaps 2352(pA), rA3
            #endif
         #endif
         #if KB > 38
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 39
               movddup 184(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 39
               movaps 2368(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 39
               movaps 2384(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 39
               movaps 2400(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 39
               movaps 2416(pA), rA3
            #endif
         #endif
         #if KB > 39
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 40
               movddup 192(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 40
               movaps 2432(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 40
               movaps 2448(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 40
               movaps 2464(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 40
               movaps 2480(pA), rA3
            #endif
         #endif
         #if KB > 40
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 41
               movddup 200(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 41
               movaps 2496(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 41
               movaps 2512(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 41
               movaps 2528(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 41
               movaps 2544(pA), rA3
            #endif
         #endif
         #if KB > 41
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 42
               movddup 208(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 42
               movaps 2560(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 42
               movaps 2576(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 42
               movaps 2592(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 42
               movaps 2608(pA), rA3
            #endif
         #endif
         #if KB > 42
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 43
               movddup 216(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 43
               movaps 2624(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 43
               movaps 2640(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 43
               movaps 2656(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 43
               movaps 2672(pA), rA3
            #endif
         #endif
         #if KB > 43
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 44
               movddup 224(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 44
               movaps 2688(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 44
               movaps 2704(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 44
               movaps 2720(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 44
               movaps 2736(pA), rA3
            #endif
         #endif
         #if KB > 44
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 45
               movddup 232(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 45
               movaps 2752(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 45
               movaps 2768(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 45
               movaps 2784(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 45
               movaps 2800(pA), rA3
            #endif
         #endif
         #if KB > 45
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 46
               movddup 240(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 46
               movaps 2816(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 46
               movaps 2832(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 46
               movaps 2848(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 46
               movaps 2864(pA), rA3
            #endif
         #endif
         #if KB > 46
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 47
               movddup 248(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 47
               movaps 2880(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 47
               movaps 2896(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 47
               movaps 2912(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 47
               movaps 2928(pA), rA3
            #endif
         #endif
         #if KB > 47
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 48
               movddup 256(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 48
               movaps 2944(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 48
               movaps 2960(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 48
               movaps 2976(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 48
               movaps 2992(pA), rA3
            #endif
         #endif
         #if KB > 48
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 49
               movddup 264(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 49
               movaps 3008(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 49
               movaps 3024(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 49
               movaps 3040(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 49
               movaps 3056(pA), rA3
            #endif
         #endif
         #if KB > 49
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 50
               movddup 272(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 50
               movaps 3072(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 50
               movaps 3088(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 50
               movaps 3104(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 50
               movaps 3120(pA), rA3
            #endif
         #endif
         #if KB > 50
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 51
               movddup 280(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 51
               movaps 3136(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 51
               movaps 3152(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 51
               movaps 3168(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 51
               movaps 3184(pA), rA3
            #endif
         #endif
         #if KB > 51
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 52
               movddup 288(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 52
               movaps 3200(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 52
               movaps 3216(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 52
               movaps 3232(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 52
               movaps 3248(pA), rA3
            #endif
         #endif
         #if KB > 52
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 53
               movddup 296(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 53
               movaps 3264(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 53
               movaps 3280(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 53
               movaps 3296(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 53
               movaps 3312(pA), rA3
            #endif
         #endif
         #if KB > 53
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 54
               movddup 304(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 54
               movaps 3328(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 54
               movaps 3344(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 54
               movaps 3360(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 54
               movaps 3376(pA), rA3
            #endif
         #endif
         #if KB > 54
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 55
               movddup 312(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 55
               movaps 3392(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 55
               movaps 3408(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 55
               movaps 3424(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 55
               movaps 3440(pA), rA3
            #endif
         #endif
         #if KB > 55
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 56
               movddup 320(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 56
               movaps 3456(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 56
               movaps 3472(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 56
               movaps 3488(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 56
               movaps 3504(pA), rA3
            #endif
         #endif
         #if KB > 56
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 57
               movddup 328(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 57
               movaps 3520(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 57
               movaps 3536(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 57
               movaps 3552(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 57
               movaps 3568(pA), rA3
            #endif
         #endif
         #if KB > 57
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 58
               movddup 336(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 58
               movaps 3584(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 58
               movaps 3600(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 58
               movaps 3616(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 58
               movaps 3632(pA), rA3
            #endif
         #endif
         #if KB > 58
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 59
               movddup 344(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 59
               movaps 3648(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 59
               movaps 3664(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 59
               movaps 3680(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 59
               movaps 3696(pA), rA3
            #endif
         #endif
         #if KB > 59
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 60
               movddup 352(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 60
               movaps 3712(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 60
               movaps 3728(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 60
               movaps 3744(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 60
               movaps 3760(pA), rA3
            #endif
         #endif
         #if KB > 60
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 61
               movddup 360(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 61
               movaps 3776(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 61
               movaps 3792(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 61
               movaps 3808(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 61
               movaps 3824(pA), rA3
            #endif
         #endif
         #if KB > 61
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 62
               movddup 368(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 62
               movaps 3840(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 62
               movaps 3856(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 62
               movaps 3872(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 62
               movaps 3888(pA), rA3
            #endif
         #endif
         #if KB > 62
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 63
               movddup 376(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 63
               movaps 3904(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 63
               movaps 3920(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 63
               movaps 3936(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 63
               movaps 3952(pA), rA3
            #endif
         #endif
         #if KB > 63
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 64
               movddup 384(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 64
               movaps 3968(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 64
               movaps 3984(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 64
               movaps 4000(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 64
               movaps 4016(pA), rA3
            #endif
         #endif
         #if KB > 64
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 65
               movddup 392(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 65
               movaps 4032(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 65
               movaps 4048(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 65
               movaps 4064(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 65
               movaps 4080(pA), rA3
            #endif
         #endif
         #if KB > 65
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 66
               movddup 400(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 66
               movaps 4096(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 66
               movaps 4112(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 66
               movaps 4128(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 66
               movaps 4144(pA), rA3
            #endif
         #endif
         #if KB > 66
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 67
               movddup 408(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 67
               movaps 4160(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 67
               movaps 4176(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 67
               movaps 4192(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 67
               movaps 4208(pA), rA3
            #endif
         #endif
         #if KB > 67
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 68
               movddup 416(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 68
               movaps 4224(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 68
               movaps 4240(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 68
               movaps 4256(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 68
               movaps 4272(pA), rA3
            #endif
         #endif
         #if KB > 68
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 69
               movddup 424(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 69
               movaps 4288(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 69
               movaps 4304(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 69
               movaps 4320(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 69
               movaps 4336(pA), rA3
            #endif
         #endif
         #if KB > 69
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 70
               movddup 432(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 70
               movaps 4352(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 70
               movaps 4368(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 70
               movaps 4384(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 70
               movaps 4400(pA), rA3
            #endif
         #endif
         #if KB > 70
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 71
               movddup 440(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 71
               movaps 4416(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 71
               movaps 4432(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 71
               movaps 4448(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 71
               movaps 4464(pA), rA3
            #endif
         #endif
         #if KB > 71
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 72
               movddup 448(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 72
               movaps 4480(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 72
               movaps 4496(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 72
               movaps 4512(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 72
               movaps 4528(pA), rA3
            #endif
         #endif
         #if KB > 72
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 73
               movddup 456(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 73
               movaps 4544(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 73
               movaps 4560(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 73
               movaps 4576(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 73
               movaps 4592(pA), rA3
            #endif
         #endif
         #if KB > 73
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 74
               movddup 464(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 74
               movaps 4608(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 74
               movaps 4624(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 74
               movaps 4640(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 74
               movaps 4656(pA), rA3
            #endif
         #endif
         #if KB > 74
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 75
               movddup 472(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 75
               movaps 4672(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 75
               movaps 4688(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 75
               movaps 4704(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 75
               movaps 4720(pA), rA3
            #endif
         #endif
         #if KB > 75
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 76
               movddup 480(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 76
               movaps 4736(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 76
               movaps 4752(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 76
               movaps 4768(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 76
               movaps 4784(pA), rA3
            #endif
         #endif
         #if KB > 76
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 77
               movddup 488(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 77
               movaps 4800(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 77
               movaps 4816(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 77
               movaps 4832(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 77
               movaps 4848(pA), rA3
            #endif
         #endif
         #if KB > 77
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 78
               movddup 496(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 78
               movaps 4864(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 78
               movaps 4880(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 78
               movaps 4896(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 78
               movaps 4912(pA), rA3
            #endif
         #endif
         #if KB > 78
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 79
               movddup 504(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 79
               movaps 4928(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 79
               movaps 4944(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 79
               movaps 4960(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 79
               movaps 4976(pA), rA3
            #endif
         #endif
         #if KB > 79
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 80
               movddup 512(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 80
               movaps 4992(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 80
               movaps 5008(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 80
               movaps 5024(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 80
               movaps 5040(pA), rA3
            #endif
         #endif
         #if KB > 80
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 81
               movddup 520(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 81
               movaps 5056(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 81
               movaps 5072(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 81
               movaps 5088(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 81
               movaps 5104(pA), rA3
            #endif
         #endif
         #if KB > 81
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 82
               movddup 528(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 82
               movaps 5120(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 82
               movaps 5136(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 82
               movaps 5152(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 82
               movaps 5168(pA), rA3
            #endif
         #endif
         #if KB > 82
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 83
               movddup 536(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 83
               movaps 5184(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 83
               movaps 5200(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 83
               movaps 5216(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 83
               movaps 5232(pA), rA3
            #endif
         #endif
         #if KB > 83
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 84
               movddup 544(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 84
               movaps 5248(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 84
               movaps 5264(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 84
               movaps 5280(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 84
               movaps 5296(pA), rA3
            #endif
         #endif
         #if KB > 84
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 85
               movddup 552(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 85
               movaps 5312(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 85
               movaps 5328(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 85
               movaps 5344(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 85
               movaps 5360(pA), rA3
            #endif
         #endif
         #if KB > 85
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 86
               movddup 560(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 86
               movaps 5376(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 86
               movaps 5392(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 86
               movaps 5408(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 86
               movaps 5424(pA), rA3
            #endif
         #endif
         #if KB > 86
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 87
               movddup 568(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 87
               movaps 5440(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 87
               movaps 5456(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 87
               movaps 5472(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 87
               movaps 5488(pA), rA3
            #endif
         #endif
         #if KB > 87
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 88
               movddup 576(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 88
               movaps 5504(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 88
               movaps 5520(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 88
               movaps 5536(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 88
               movaps 5552(pA), rA3
            #endif
         #endif
         #if KB > 88
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 89
               movddup 584(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 89
               movaps 5568(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 89
               movaps 5584(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 89
               movaps 5600(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 89
               movaps 5616(pA), rA3
            #endif
         #endif
         #if KB > 89
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 90
               movddup 592(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 90
               movaps 5632(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 90
               movaps 5648(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 90
               movaps 5664(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 90
               movaps 5680(pA), rA3
            #endif
         #endif
         #if KB > 90
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 91
               movddup 600(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 91
               movaps 5696(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 91
               movaps 5712(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 91
               movaps 5728(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 91
               movaps 5744(pA), rA3
            #endif
         #endif
         #if KB > 91
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 92
               movddup 608(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 92
               movaps 5760(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 92
               movaps 5776(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 92
               movaps 5792(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 92
               movaps 5808(pA), rA3
            #endif
         #endif
         #if KB > 92
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 93
               movddup 616(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 93
               movaps 5824(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 93
               movaps 5840(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 93
               movaps 5856(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 93
               movaps 5872(pA), rA3
            #endif
         #endif
         #if KB > 93
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 94
               movddup 624(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 94
               movaps 5888(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 94
               movaps 5904(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 94
               movaps 5920(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 94
               movaps 5936(pA), rA3
            #endif
         #endif
         #if KB > 94
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 95
               movddup 632(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 95
               movaps 5952(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 95
               movaps 5968(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 95
               movaps 5984(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 95
               movaps 6000(pA), rA3
            #endif
         #endif
         #if KB > 95
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 96
               movddup 640(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 96
               movaps 6016(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 96
               movaps 6032(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 96
               movaps 6048(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 96
               movaps 6064(pA), rA3
            #endif
         #endif
         #if KB > 96
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 97
               movddup 648(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 97
               movaps 6080(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 97
               movaps 6096(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 97
               movaps 6112(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 97
               movaps 6128(pA), rA3
            #endif
         #endif
         #if KB > 97
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 98
               movddup 656(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 98
               movaps 6144(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 98
               movaps 6160(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 98
               movaps 6176(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 98
               movaps 6192(pA), rA3
            #endif
         #endif
         #if KB > 98
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 99
               movddup 664(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 99
               movaps 6208(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 99
               movaps 6224(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 99
               movaps 6240(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 99
               movaps 6256(pA), rA3
            #endif
         #endif
         #if KB > 99
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 100
               movddup 672(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 100
               movaps 6272(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 100
               movaps 6288(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 100
               movaps 6304(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 100
               movaps 6320(pA), rA3
            #endif
         #endif
         #if KB > 100
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 101
               movddup 680(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 101
               movaps 6336(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 101
               movaps 6352(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 101
               movaps 6368(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 101
               movaps 6384(pA), rA3
            #endif
         #endif
         #if KB > 101
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 102
               movddup 688(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 102
               movaps 6400(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 102
               movaps 6416(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 102
               movaps 6432(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 102
               movaps 6448(pA), rA3
            #endif
         #endif
         #if KB > 102
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 103
               movddup 696(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 103
               movaps 6464(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 103
               movaps 6480(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 103
               movaps 6496(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 103
               movaps 6512(pA), rA3
            #endif
         #endif
         #if KB > 103
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 104
               movddup 704(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 104
               movaps 6528(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 104
               movaps 6544(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 104
               movaps 6560(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 104
               movaps 6576(pA), rA3
            #endif
         #endif
         #if KB > 104
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 105
               movddup 712(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 105
               movaps 6592(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 105
               movaps 6608(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 105
               movaps 6624(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 105
               movaps 6640(pA), rA3
            #endif
         #endif
         #if KB > 105
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 106
               movddup 720(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 106
               movaps 6656(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 106
               movaps 6672(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 106
               movaps 6688(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 106
               movaps 6704(pA), rA3
            #endif
         #endif
         #if KB > 106
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 107
               movddup 728(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 107
               movaps 6720(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 107
               movaps 6736(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 107
               movaps 6752(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 107
               movaps 6768(pA), rA3
            #endif
         #endif
         #if KB > 107
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 108
               movddup 736(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 108
               movaps 6784(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 108
               movaps 6800(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 108
               movaps 6816(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 108
               movaps 6832(pA), rA3
            #endif
         #endif
         #if KB > 108
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 109
               movddup 744(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 109
               movaps 6848(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 109
               movaps 6864(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 109
               movaps 6880(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 109
               movaps 6896(pA), rA3
            #endif
         #endif
         #if KB > 109
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 110
               movddup 752(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 110
               movaps 6912(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 110
               movaps 6928(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 110
               movaps 6944(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 110
               movaps 6960(pA), rA3
            #endif
         #endif
         #if KB > 110
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 111
               movddup 760(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 111
               movaps 6976(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 111
               movaps 6992(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 111
               movaps 7008(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 111
               movaps 7024(pA), rA3
            #endif
         #endif
         #if KB > 111
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 112
               movddup 768(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 112
               movaps 7040(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 112
               movaps 7056(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 112
               movaps 7072(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 112
               movaps 7088(pA), rA3
            #endif
         #endif
         #if KB > 112
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 113
               movddup 776(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 113
               movaps 7104(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 113
               movaps 7120(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 113
               movaps 7136(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 113
               movaps 7152(pA), rA3
            #endif
         #endif
         #if KB > 113
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 114
               movddup 784(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 114
               movaps 7168(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 114
               movaps 7184(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 114
               movaps 7200(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 114
               movaps 7216(pA), rA3
            #endif
         #endif
         #if KB > 114
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 115
               movddup 792(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 115
               movaps 7232(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 115
               movaps 7248(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 115
               movaps 7264(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 115
               movaps 7280(pA), rA3
            #endif
         #endif
         #if KB > 115
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 116
               movddup 800(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 116
               movaps 7296(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 116
               movaps 7312(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 116
               movaps 7328(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 116
               movaps 7344(pA), rA3
            #endif
         #endif
         #if KB > 116
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 117
               movddup 808(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 117
               movaps 7360(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 117
               movaps 7376(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 117
               movaps 7392(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 117
               movaps 7408(pA), rA3
            #endif
         #endif
         #if KB > 117
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 118
               movddup 816(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 118
               movaps 7424(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 118
               movaps 7440(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 118
               movaps 7456(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 118
               movaps 7472(pA), rA3
            #endif
         #endif
         #if KB > 118
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 119
               movddup 824(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 119
               movaps 7488(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 119
               movaps 7504(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 119
               movaps 7520(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 119
               movaps 7536(pA), rA3
            #endif
         #endif
         #if KB > 119
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 120
               movddup 832(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 120
               movaps 7552(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 120
               movaps 7568(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 120
               movaps 7584(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 120
               movaps 7600(pA), rA3
            #endif
         #endif
         #if KB > 120
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 121
               movddup 840(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 121
               movaps 7616(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 121
               movaps 7632(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 121
               movaps 7648(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 121
               movaps 7664(pA), rA3
            #endif
         #endif
         #if KB > 121
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 122
               movddup 848(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 122
               movaps 7680(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 122
               movaps 7696(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 122
               movaps 7712(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 122
               movaps 7728(pA), rA3
            #endif
         #endif
         #if KB > 122
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 123
               movddup 856(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 123
               movaps 7744(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 123
               movaps 7760(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 123
               movaps 7776(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 123
               movaps 7792(pA), rA3
            #endif
         #endif
         #if KB > 123
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 124
               movddup 864(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 124
               movaps 7808(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 124
               movaps 7824(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 124
               movaps 7840(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 124
               movaps 7856(pA), rA3
            #endif
         #endif
         #if KB > 124
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 125
               movddup 872(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 125
               movaps 7872(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 125
               movaps 7888(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 125
               movaps 7904(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 125
               movaps 7920(pA), rA3
            #endif
         #endif
         #if KB > 125
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 126
               movddup 880(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 126
               movaps 7936(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 126
               movaps 7952(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 126
               movaps 7968(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 126
               movaps 7984(pA), rA3
            #endif
         #endif
         #if KB > 126
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 127
               movddup 888(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 127
               movaps 8000(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 127
               movaps 8016(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 127
               movaps 8032(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 127
               movaps 8048(pA), rA3
            #endif
         #endif
         #if KB > 127
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 128
               movddup 896(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 128
               movaps 8064(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 128
               movaps 8080(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 128
               movaps 8096(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 128
               movaps 8112(pA), rA3
            #endif
         #endif
         #if KB > 128
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 129
               movddup 904(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 129
               movaps 8128(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 129
               movaps 8144(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 129
               movaps 8160(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 129
               movaps 8176(pA), rA3
            #endif
         #endif
         #if KB > 129
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 130
               movddup 912(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 130
               movaps 8192(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 130
               movaps 8208(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 130
               movaps 8224(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 130
               movaps 8240(pA), rA3
            #endif
         #endif
         #if KB > 130
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 131
               movddup 920(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 131
               movaps 8256(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 131
               movaps 8272(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 131
               movaps 8288(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 131
               movaps 8304(pA), rA3
            #endif
         #endif
         #if KB > 131
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 132
               movddup 928(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 132
               movaps 8320(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 132
               movaps 8336(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 132
               movaps 8352(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 132
               movaps 8368(pA), rA3
            #endif
         #endif
         #if KB > 132
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 133
               movddup 936(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 133
               movaps 8384(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 133
               movaps 8400(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 133
               movaps 8416(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 133
               movaps 8432(pA), rA3
            #endif
         #endif
         #if KB > 133
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 134
               movddup 944(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 134
               movaps 8448(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 134
               movaps 8464(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 134
               movaps 8480(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 134
               movaps 8496(pA), rA3
            #endif
         #endif
         #if KB > 134
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 135
               movddup 952(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 135
               movaps 8512(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 135
               movaps 8528(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 135
               movaps 8544(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 135
               movaps 8560(pA), rA3
            #endif
         #endif
         #if KB > 135
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 136
               movddup 960(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 136
               movaps 8576(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 136
               movaps 8592(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 136
               movaps 8608(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 136
               movaps 8624(pA), rA3
            #endif
         #endif
         #if KB > 136
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 137
               movddup 968(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 137
               movaps 8640(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 137
               movaps 8656(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 137
               movaps 8672(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 137
               movaps 8688(pA), rA3
            #endif
         #endif
         #if KB > 137
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 138
               movddup 976(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 138
               movaps 8704(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 138
               movaps 8720(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 138
               movaps 8736(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 138
               movaps 8752(pA), rA3
            #endif
         #endif
         #if KB > 138
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 139
               movddup 984(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 139
               movaps 8768(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 139
               movaps 8784(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 139
               movaps 8800(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 139
               movaps 8816(pA), rA3
            #endif
         #endif
         #if KB > 139
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 140
               movddup 992(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 140
               movaps 8832(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 140
               movaps 8848(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 140
               movaps 8864(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 140
               movaps 8880(pA), rA3
            #endif
         #endif
         #if KB > 140
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 141
               movddup 1000(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 141
               movaps 8896(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 141
               movaps 8912(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 141
               movaps 8928(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 141
               movaps 8944(pA), rA3
            #endif
         #endif
         #if KB > 141
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 142
               movddup 1008(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 142
               movaps 8960(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 142
               movaps 8976(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 142
               movaps 8992(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 142
               movaps 9008(pA), rA3
            #endif
         #endif
         #if KB > 142
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 143
               movddup 1016(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 143
               movaps 9024(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 143
               movaps 9040(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 143
               movaps 9056(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 143
               movaps 9072(pA), rA3
            #endif
         #endif
         #if KB > 143
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 144
               movddup 1024(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 144
               movaps 9088(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 144
               movaps 9104(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 144
               movaps 9120(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 144
               movaps 9136(pA), rA3
            #endif
         #endif
         #if KB > 144
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 145
               movddup 1032(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 145
               movaps 9152(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 145
               movaps 9168(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 145
               movaps 9184(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 145
               movaps 9200(pA), rA3
            #endif
         #endif
         #if KB > 145
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 146
               movddup 1040(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 146
               movaps 9216(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 146
               movaps 9232(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 146
               movaps 9248(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 146
               movaps 9264(pA), rA3
            #endif
         #endif
         #if KB > 146
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 147
               movddup 1048(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 147
               movaps 9280(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 147
               movaps 9296(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 147
               movaps 9312(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 147
               movaps 9328(pA), rA3
            #endif
         #endif
         #if KB > 147
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 148
               movddup 1056(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 148
               movaps 9344(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 148
               movaps 9360(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 148
               movaps 9376(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 148
               movaps 9392(pA), rA3
            #endif
         #endif
         #if KB > 148
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 149
               movddup 1064(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 149
               movaps 9408(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 149
               movaps 9424(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 149
               movaps 9440(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 149
               movaps 9456(pA), rA3
            #endif
         #endif
         #if KB > 149
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 150
               movddup 1072(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 150
               movaps 9472(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 150
               movaps 9488(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 150
               movaps 9504(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 150
               movaps 9520(pA), rA3
            #endif
         #endif
         #if KB > 150
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 151
               movddup 1080(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 151
               movaps 9536(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 151
               movaps 9552(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 151
               movaps 9568(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 151
               movaps 9584(pA), rA3
            #endif
         #endif
         #if KB > 151
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 152
               movddup 1088(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 152
               movaps 9600(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 152
               movaps 9616(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 152
               movaps 9632(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 152
               movaps 9648(pA), rA3
            #endif
         #endif
         #if KB > 152
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 153
               movddup 1096(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 153
               movaps 9664(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 153
               movaps 9680(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 153
               movaps 9696(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 153
               movaps 9712(pA), rA3
            #endif
         #endif
         #if KB > 153
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 154
               movddup 1104(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 154
               movaps 9728(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 154
               movaps 9744(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 154
               movaps 9760(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 154
               movaps 9776(pA), rA3
            #endif
         #endif
         #if KB > 154
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 155
               movddup 1112(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 155
               movaps 9792(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 155
               movaps 9808(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 155
               movaps 9824(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 155
               movaps 9840(pA), rA3
            #endif
         #endif
         #if KB > 155
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 156
               movddup 1120(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 156
               movaps 9856(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 156
               movaps 9872(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 156
               movaps 9888(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 156
               movaps 9904(pA), rA3
            #endif
         #endif
         #if KB > 156
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 157
               movddup 1128(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 157
               movaps 9920(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 157
               movaps 9936(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 157
               movaps 9952(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 157
               movaps 9968(pA), rA3
            #endif
         #endif
         #if KB > 157
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 158
               movddup 1136(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 158
               movaps 9984(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 158
               movaps 10000(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 158
               movaps 10016(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 158
               movaps 10032(pA), rA3
            #endif
         #endif
         #if KB > 158
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 159
               movddup 1144(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 159
               movaps 10048(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 159
               movaps 10064(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 159
               movaps 10080(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 159
               movaps 10096(pA), rA3
            #endif
         #endif
         #if KB > 159
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 160
               movddup 1152(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 160
               movaps 10112(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 160
               movaps 10128(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 160
               movaps 10144(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 160
               movaps 10160(pA), rA3
            #endif
         #endif
         #if KB > 160
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 161
               movddup 1160(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 161
               movaps 10176(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 161
               movaps 10192(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 161
               movaps 10208(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 161
               movaps 10224(pA), rA3
            #endif
         #endif
         #if KB > 161
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 162
               movddup 1168(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 162
               movaps 10240(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 162
               movaps 10256(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 162
               movaps 10272(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 162
               movaps 10288(pA), rA3
            #endif
         #endif
         #if KB > 162
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 163
               movddup 1176(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 163
               movaps 10304(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 163
               movaps 10320(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 163
               movaps 10336(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 163
               movaps 10352(pA), rA3
            #endif
         #endif
         #if KB > 163
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 164
               movddup 1184(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 164
               movaps 10368(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 164
               movaps 10384(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 164
               movaps 10400(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 164
               movaps 10416(pA), rA3
            #endif
         #endif
         #if KB > 164
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 165
               movddup 1192(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 165
               movaps 10432(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 165
               movaps 10448(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 165
               movaps 10464(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 165
               movaps 10480(pA), rA3
            #endif
         #endif
         #if KB > 165
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 166
               movddup 1200(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 166
               movaps 10496(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 166
               movaps 10512(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 166
               movaps 10528(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 166
               movaps 10544(pA), rA3
            #endif
         #endif
         #if KB > 166
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 167
               movddup 1208(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 167
               movaps 10560(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 167
               movaps 10576(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 167
               movaps 10592(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 167
               movaps 10608(pA), rA3
            #endif
         #endif
         #if KB > 167
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 168
               movddup 1216(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 168
               movaps 10624(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 168
               movaps 10640(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 168
               movaps 10656(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 168
               movaps 10672(pA), rA3
            #endif
         #endif
         #if KB > 168
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 169
               movddup 1224(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 169
               movaps 10688(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 169
               movaps 10704(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 169
               movaps 10720(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 169
               movaps 10736(pA), rA3
            #endif
         #endif
         #if KB > 169
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 170
               movddup 1232(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 170
               movaps 10752(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 170
               movaps 10768(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 170
               movaps 10784(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 170
               movaps 10800(pA), rA3
            #endif
         #endif
         #if KB > 170
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 171
               movddup 1240(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 171
               movaps 10816(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 171
               movaps 10832(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 171
               movaps 10848(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 171
               movaps 10864(pA), rA3
            #endif
         #endif
         #if KB > 171
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 172
               movddup 1248(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 172
               movaps 10880(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 172
               movaps 10896(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 172
               movaps 10912(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 172
               movaps 10928(pA), rA3
            #endif
         #endif
         #if KB > 172
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 173
               movddup 1256(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 173
               movaps 10944(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 173
               movaps 10960(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 173
               movaps 10976(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 173
               movaps 10992(pA), rA3
            #endif
         #endif
         #if KB > 173
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 174
               movddup 1264(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 174
               movaps 11008(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 174
               movaps 11024(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 174
               movaps 11040(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 174
               movaps 11056(pA), rA3
            #endif
         #endif
         #if KB > 174
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 175
               movddup 1272(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 175
               movaps 11072(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 175
               movaps 11088(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 175
               movaps 11104(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 175
               movaps 11120(pA), rA3
            #endif
         #endif
         #if KB > 175
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 176
               movddup 1280(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 176
               movaps 11136(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 176
               movaps 11152(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 176
               movaps 11168(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 176
               movaps 11184(pA), rA3
            #endif
         #endif
         #if KB > 176
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 177
               movddup 1288(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 177
               movaps 11200(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 177
               movaps 11216(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 177
               movaps 11232(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 177
               movaps 11248(pA), rA3
            #endif
         #endif
         #if KB > 177
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 178
               movddup 1296(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 178
               movaps 11264(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 178
               movaps 11280(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 178
               movaps 11296(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 178
               movaps 11312(pA), rA3
            #endif
         #endif
         #if KB > 178
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 179
               movddup 1304(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 179
               movaps 11328(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 179
               movaps 11344(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 179
               movaps 11360(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 179
               movaps 11376(pA), rA3
            #endif
         #endif
         #if KB > 179
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 180
               movddup 1312(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 180
               movaps 11392(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 180
               movaps 11408(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 180
               movaps 11424(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 180
               movaps 11440(pA), rA3
            #endif
         #endif
         #if KB > 180
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 181
               movddup 1320(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 181
               movaps 11456(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 181
               movaps 11472(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 181
               movaps 11488(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 181
               movaps 11504(pA), rA3
            #endif
         #endif
         #if KB > 181
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 182
               movddup 1328(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 182
               movaps 11520(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 182
               movaps 11536(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 182
               movaps 11552(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 182
               movaps 11568(pA), rA3
            #endif
         #endif
         #if KB > 182
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 183
               movddup 1336(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 183
               movaps 11584(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 183
               movaps 11600(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 183
               movaps 11616(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 183
               movaps 11632(pA), rA3
            #endif
         #endif
         #if KB > 183
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 184
               movddup 1344(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 184
               movaps 11648(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 184
               movaps 11664(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 184
               movaps 11680(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 184
               movaps 11696(pA), rA3
            #endif
         #endif
         #if KB > 184
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 185
               movddup 1352(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 185
               movaps 11712(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 185
               movaps 11728(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 185
               movaps 11744(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 185
               movaps 11760(pA), rA3
            #endif
         #endif
         #if KB > 185
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 186
               movddup 1360(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 186
               movaps 11776(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 186
               movaps 11792(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 186
               movaps 11808(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 186
               movaps 11824(pA), rA3
            #endif
         #endif
         #if KB > 186
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 187
               movddup 1368(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 187
               movaps 11840(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 187
               movaps 11856(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 187
               movaps 11872(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 187
               movaps 11888(pA), rA3
            #endif
         #endif
         #if KB > 187
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 188
               movddup 1376(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 188
               movaps 11904(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 188
               movaps 11920(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 188
               movaps 11936(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 188
               movaps 11952(pA), rA3
            #endif
         #endif
         #if KB > 188
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 189
               movddup 1384(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 189
               movaps 11968(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 189
               movaps 11984(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 189
               movaps 12000(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 189
               movaps 12016(pA), rA3
            #endif
         #endif
         #if KB > 189
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 190
               movddup 1392(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 190
               movaps 12032(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 190
               movaps 12048(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 190
               movaps 12064(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 190
               movaps 12080(pA), rA3
            #endif
         #endif
         #if KB > 190
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 191
               movddup 1400(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 191
               movaps 12096(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 191
               movaps 12112(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 191
               movaps 12128(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 191
               movaps 12144(pA), rA3
            #endif
         #endif
         #if KB > 191
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 192
               movddup 1408(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 192
               movaps 12160(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 192
               movaps 12176(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 192
               movaps 12192(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 192
               movaps 12208(pA), rA3
            #endif
         #endif
         #if KB > 192
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 193
               movddup 1416(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 193
               movaps 12224(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 193
               movaps 12240(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 193
               movaps 12256(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 193
               movaps 12272(pA), rA3
            #endif
         #endif
         #if KB > 193
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 194
               movddup 1424(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 194
               movaps 12288(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 194
               movaps 12304(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 194
               movaps 12320(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 194
               movaps 12336(pA), rA3
            #endif
         #endif
         #if KB > 194
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 195
               movddup 1432(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 195
               movaps 12352(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 195
               movaps 12368(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 195
               movaps 12384(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 195
               movaps 12400(pA), rA3
            #endif
         #endif
         #if KB > 195
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 196
               movddup 1440(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 196
               movaps 12416(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 196
               movaps 12432(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 196
               movaps 12448(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 196
               movaps 12464(pA), rA3
            #endif
         #endif
         #if KB > 196
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 197
               movddup 1448(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 197
               movaps 12480(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 197
               movaps 12496(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 197
               movaps 12512(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 197
               movaps 12528(pA), rA3
            #endif
         #endif
         #if KB > 197
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 198
               movddup 1456(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 198
               movaps 12544(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 198
               movaps 12560(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 198
               movaps 12576(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 198
               movaps 12592(pA), rA3
            #endif
         #endif
         #if KB > 198
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 199
               movddup 1464(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 199
               movaps 12608(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 199
               movaps 12624(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 199
               movaps 12640(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 199
               movaps 12656(pA), rA3
            #endif
         #endif
         #if KB > 199
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 200
               movddup 1472(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 200
               movaps 12672(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 200
               movaps 12688(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 200
               movaps 12704(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 200
               movaps 12720(pA), rA3
            #endif
         #endif
         #if KB > 200
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 201
               movddup 1480(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 201
               movaps 12736(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 201
               movaps 12752(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 201
               movaps 12768(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 201
               movaps 12784(pA), rA3
            #endif
         #endif
         #if KB > 201
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 202
               movddup 1488(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 202
               movaps 12800(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 202
               movaps 12816(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 202
               movaps 12832(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 202
               movaps 12848(pA), rA3
            #endif
         #endif
         #if KB > 202
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 203
               movddup 1496(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 203
               movaps 12864(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 203
               movaps 12880(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 203
               movaps 12896(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 203
               movaps 12912(pA), rA3
            #endif
         #endif
         #if KB > 203
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 204
               movddup 1504(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 204
               movaps 12928(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 204
               movaps 12944(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 204
               movaps 12960(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 204
               movaps 12976(pA), rA3
            #endif
         #endif
         #if KB > 204
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 205
               movddup 1512(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 205
               movaps 12992(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 205
               movaps 13008(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 205
               movaps 13024(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 205
               movaps 13040(pA), rA3
            #endif
         #endif
         #if KB > 205
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 206
               movddup 1520(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 206
               movaps 13056(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 206
               movaps 13072(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 206
               movaps 13088(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 206
               movaps 13104(pA), rA3
            #endif
         #endif
         #if KB > 206
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 207
               movddup 1528(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 207
               movaps 13120(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 207
               movaps 13136(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 207
               movaps 13152(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 207
               movaps 13168(pA), rA3
            #endif
         #endif
         #if KB > 207
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 208
               movddup 1536(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 208
               movaps 13184(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 208
               movaps 13200(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 208
               movaps 13216(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 208
               movaps 13232(pA), rA3
            #endif
         #endif
         #if KB > 208
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 209
               movddup 1544(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 209
               movaps 13248(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 209
               movaps 13264(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 209
               movaps 13280(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 209
               movaps 13296(pA), rA3
            #endif
         #endif
         #if KB > 209
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 210
               movddup 1552(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 210
               movaps 13312(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 210
               movaps 13328(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 210
               movaps 13344(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 210
               movaps 13360(pA), rA3
            #endif
         #endif
         #if KB > 210
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 211
               movddup 1560(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 211
               movaps 13376(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 211
               movaps 13392(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 211
               movaps 13408(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 211
               movaps 13424(pA), rA3
            #endif
         #endif
         #if KB > 211
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 212
               movddup 1568(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 212
               movaps 13440(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 212
               movaps 13456(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 212
               movaps 13472(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 212
               movaps 13488(pA), rA3
            #endif
         #endif
         #if KB > 212
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 213
               movddup 1576(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 213
               movaps 13504(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 213
               movaps 13520(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 213
               movaps 13536(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 213
               movaps 13552(pA), rA3
            #endif
         #endif
         #if KB > 213
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 214
               movddup 1584(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 214
               movaps 13568(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 214
               movaps 13584(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 214
               movaps 13600(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 214
               movaps 13616(pA), rA3
            #endif
         #endif
         #if KB > 214
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 215
               movddup 1592(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 215
               movaps 13632(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 215
               movaps 13648(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 215
               movaps 13664(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 215
               movaps 13680(pA), rA3
            #endif
         #endif
         #if KB > 215
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 216
               movddup 1600(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 216
               movaps 13696(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 216
               movaps 13712(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 216
               movaps 13728(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 216
               movaps 13744(pA), rA3
            #endif
         #endif
         #if KB > 216
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 217
               movddup 1608(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 217
               movaps 13760(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 217
               movaps 13776(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 217
               movaps 13792(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 217
               movaps 13808(pA), rA3
            #endif
         #endif
         #if KB > 217
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 218
               movddup 1616(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 218
               movaps 13824(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 218
               movaps 13840(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 218
               movaps 13856(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 218
               movaps 13872(pA), rA3
            #endif
         #endif
         #if KB > 218
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 219
               movddup 1624(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 219
               movaps 13888(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 219
               movaps 13904(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 219
               movaps 13920(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 219
               movaps 13936(pA), rA3
            #endif
         #endif
         #if KB > 219
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 220
               movddup 1632(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 220
               movaps 13952(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 220
               movaps 13968(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 220
               movaps 13984(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 220
               movaps 14000(pA), rA3
            #endif
         #endif
         #if KB > 220
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 221
               movddup 1640(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 221
               movaps 14016(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 221
               movaps 14032(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 221
               movaps 14048(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 221
               movaps 14064(pA), rA3
            #endif
         #endif
         #if KB > 221
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 222
               movddup 1648(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 222
               movaps 14080(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 222
               movaps 14096(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 222
               movaps 14112(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 222
               movaps 14128(pA), rA3
            #endif
         #endif
         #if KB > 222
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 223
               movddup 1656(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 223
               movaps 14144(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 223
               movaps 14160(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 223
               movaps 14176(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 223
               movaps 14192(pA), rA3
            #endif
         #endif
         #if KB > 223
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 224
               movddup 1664(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 224
               movaps 14208(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 224
               movaps 14224(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 224
               movaps 14240(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 224
               movaps 14256(pA), rA3
            #endif
         #endif
         #if KB > 224
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 225
               movddup 1672(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 225
               movaps 14272(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 225
               movaps 14288(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 225
               movaps 14304(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 225
               movaps 14320(pA), rA3
            #endif
         #endif
         #if KB > 225
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 226
               movddup 1680(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 226
               movaps 14336(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 226
               movaps 14352(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 226
               movaps 14368(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 226
               movaps 14384(pA), rA3
            #endif
         #endif
         #if KB > 226
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 227
               movddup 1688(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 227
               movaps 14400(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 227
               movaps 14416(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 227
               movaps 14432(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 227
               movaps 14448(pA), rA3
            #endif
         #endif
         #if KB > 227
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 228
               movddup 1696(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 228
               movaps 14464(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 228
               movaps 14480(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 228
               movaps 14496(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 228
               movaps 14512(pA), rA3
            #endif
         #endif
         #if KB > 228
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 229
               movddup 1704(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 229
               movaps 14528(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 229
               movaps 14544(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 229
               movaps 14560(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 229
               movaps 14576(pA), rA3
            #endif
         #endif
         #if KB > 229
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 230
               movddup 1712(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 230
               movaps 14592(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 230
               movaps 14608(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 230
               movaps 14624(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 230
               movaps 14640(pA), rA3
            #endif
         #endif
         #if KB > 230
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 231
               movddup 1720(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 231
               movaps 14656(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 231
               movaps 14672(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 231
               movaps 14688(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 231
               movaps 14704(pA), rA3
            #endif
         #endif
         #if KB > 231
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 232
               movddup 1728(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 232
               movaps 14720(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 232
               movaps 14736(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 232
               movaps 14752(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 232
               movaps 14768(pA), rA3
            #endif
         #endif
         #if KB > 232
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 233
               movddup 1736(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 233
               movaps 14784(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 233
               movaps 14800(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 233
               movaps 14816(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 233
               movaps 14832(pA), rA3
            #endif
         #endif
         #if KB > 233
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 234
               movddup 1744(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 234
               movaps 14848(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 234
               movaps 14864(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 234
               movaps 14880(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 234
               movaps 14896(pA), rA3
            #endif
         #endif
         #if KB > 234
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 235
               movddup 1752(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 235
               movaps 14912(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 235
               movaps 14928(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 235
               movaps 14944(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 235
               movaps 14960(pA), rA3
            #endif
         #endif
         #if KB > 235
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 236
               movddup 1760(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 236
               movaps 14976(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 236
               movaps 14992(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 236
               movaps 15008(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 236
               movaps 15024(pA), rA3
            #endif
         #endif
         #if KB > 236
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 237
               movddup 1768(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 237
               movaps 15040(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 237
               movaps 15056(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 237
               movaps 15072(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 237
               movaps 15088(pA), rA3
            #endif
         #endif
         #if KB > 237
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 238
               movddup 1776(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 238
               movaps 15104(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 238
               movaps 15120(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 238
               movaps 15136(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 238
               movaps 15152(pA), rA3
            #endif
         #endif
         #if KB > 238
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 239
               movddup 1784(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 239
               movaps 15168(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 239
               movaps 15184(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 239
               movaps 15200(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 239
               movaps 15216(pA), rA3
            #endif
         #endif
         #if KB > 239
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 240
               movddup 1792(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 240
               movaps 15232(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 240
               movaps 15248(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 240
               movaps 15264(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 240
               movaps 15280(pA), rA3
            #endif
         #endif
         #if KB > 240
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 241
               movddup 1800(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 241
               movaps 15296(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 241
               movaps 15312(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 241
               movaps 15328(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 241
               movaps 15344(pA), rA3
            #endif
         #endif
         #if KB > 241
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 242
               movddup 1808(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 242
               movaps 15360(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 242
               movaps 15376(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 242
               movaps 15392(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 242
               movaps 15408(pA), rA3
            #endif
         #endif
         #if KB > 242
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 243
               movddup 1816(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 243
               movaps 15424(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 243
               movaps 15440(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 243
               movaps 15456(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 243
               movaps 15472(pA), rA3
            #endif
         #endif
         #if KB > 243
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 244
               movddup 1824(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 244
               movaps 15488(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 244
               movaps 15504(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 244
               movaps 15520(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 244
               movaps 15536(pA), rA3
            #endif
         #endif
         #if KB > 244
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 245
               movddup 1832(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 245
               movaps 15552(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 245
               movaps 15568(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 245
               movaps 15584(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 245
               movaps 15600(pA), rA3
            #endif
         #endif
         #if KB > 245
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 246
               movddup 1840(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 246
               movaps 15616(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 246
               movaps 15632(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 246
               movaps 15648(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 246
               movaps 15664(pA), rA3
            #endif
         #endif
         #if KB > 246
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 247
               movddup 1848(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 247
               movaps 15680(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 247
               movaps 15696(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 247
               movaps 15712(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 247
               movaps 15728(pA), rA3
            #endif
         #endif
         #if KB > 247
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 248
               movddup 1856(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 248
               movaps 15744(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 248
               movaps 15760(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 248
               movaps 15776(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 248
               movaps 15792(pA), rA3
            #endif
         #endif
         #if KB > 248
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 249
               movddup 1864(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 249
               movaps 15808(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 249
               movaps 15824(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 249
               movaps 15840(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 249
               movaps 15856(pA), rA3
            #endif
         #endif
         #if KB > 249
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 250
               movddup 1872(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 250
               movaps 15872(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 250
               movaps 15888(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 250
               movaps 15904(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 250
               movaps 15920(pA), rA3
            #endif
         #endif
         #if KB > 250
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 251
               movddup 1880(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 251
               movaps 15936(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 251
               movaps 15952(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 251
               movaps 15968(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 251
               movaps 15984(pA), rA3
            #endif
         #endif
         #if KB > 251
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 252
               movddup 1888(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 252
               movaps 16000(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 252
               movaps 16016(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 252
               movaps 16032(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 252
               movaps 16048(pA), rA3
            #endif
         #endif
         #if KB > 252
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 253
               movddup 1896(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 253
               movaps 16064(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 253
               movaps 16080(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 253
               movaps 16096(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 253
               movaps 16112(pA), rA3
            #endif
         #endif
         #if KB > 253
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 254
               movddup 1904(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 254
               movaps 16128(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 254
               movaps 16144(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 254
               movaps 16160(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 254
               movaps 16176(pA), rA3
            #endif
         #endif
         #if KB > 254
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 255
               movddup 1912(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 255
               movaps 16192(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 255
               movaps 16208(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 255
               movaps 16224(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 255
               movaps 16240(pA), rA3
            #endif
         #endif
         #if KB > 255
            movshdup rB0, rB1
            movsldup rB0, rB0

            movaps rB0, rm0
            mulps rA0, rm0
            addps rm0, rC00

            movaps rB0, rm0
            mulps rA1, rm0
            addps rm0, rC10

            movaps rB0, rm0
            mulps rA2, rm0
            addps rm0, rC20

            mulps rA3, rB0
            addps rB0, rC30
            #if KB > 256
               movddup 1920(pB), rB0   /* b1 b0 b1 b0 */
            #endif

            mulps rB1, rA0
            addps rA0, rC01
            #if KB > 256
               movaps 16256(pA), rA0
            #endif

            mulps rB1, rA1
            addps rA1, rC11
            #if KB > 256
               movaps 16272(pA), rA1
            #endif

            mulps rB1, rA2
            addps rA2, rC21
            #if KB > 256
               movaps 16288(pA), rA2
            #endif

            mulps rB1, rA3
            addps rA3, rC31
            #if KB > 256
               movaps 16304(pA), rA3
            #endif
         #endif

         #ifndef BETA0
            VCOP (pC), rC00
         #endif
         movaps rC00, (pC)
         #ifndef BETA0
            VCOP 16(pC), rC10
         #endif
         movaps rC10, 16(pC)
         #ifndef BETA0
            VCOP 32(pC), rC20
         #endif
         movaps rC20, 32(pC)
         #ifndef BETA0
            VCOP 48(pC), rC30
         #endif
         movaps rC30, 48(pC)
         #ifndef BETA0
            VCOP 64(pC), rC01
         #endif
         movaps rC01, 64(pC)
         #ifndef BETA0
            VCOP 80(pC), rC11
         #endif
         movaps rC11, 80(pC)
         #ifndef BETA0
            VCOP 96(pC), rC21
         #endif
         movaps rC21, 96(pC)
         #ifndef BETA0
            VCOP 112(pC), rC31
         #endif
         movaps rC31, 112(pC)

         add $KB*2*4, pB
         sub $-128, pC
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
