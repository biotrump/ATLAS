#include "atlas_asm.h"

#define rm0     %xmm0
#define rB0     %xmm1
#define rB1     %xmm2
#define rB2     %xmm3
#define rA0     %xmm4
#define rA1     %xmm5
#define rA2     %xmm6
#define rC00    %xmm7
#define rC10    %xmm8
#define rC20    %xmm9
#define rC01    %xmm10
#define rC11    %xmm11
#define rC21    %xmm12
#define rC02    %xmm13
#define rC12    %xmm14
#define rC22    %xmm15
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
// #define PFADIST 128
#define PFBDIST 1088
#define prefA(m_) prefetcht2 m_
#define prefB(m_) prefetcht2 m_
//#define prefC(m_) prefetcht0 m_
#define prefC(m_) prefetchw m_
#define FMAC vfmadd231pd   /* FMAC m256/r256, rs1, rd */
#if defined(BETAN) || defined(BETAn)
   #define BETAN1
#endif
#ifdef BETAN1
   #define VCOP subps
#else
   #define VCOP addps
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
   movq 8(%rsp), pfA       /* pfA = pAn */
   movq 16(%rsp), pfB      /* pfB = pBn */
   cmp pfB, pB
   CMOVE pfA, pfB
   mov KK, KK0
   sub $-128, pC
   sub $-128, pA
   sub $-128, pB
   mov $24, r24
   mov pA, pA0
   mov pB, pB0
   mov $6*3*8, incPF
/*
 * incAm = 12*sizeof*K = 12*4*K = 16*3*K
 */
   lea (KK, KK,2), incAm   /* incAm = 3*K */
   shl $4, incAm           /* incAm = 16*3*K */

   ALIGN16
   MLOOP:
      NLOOP:
/*
 *       First iteration peeled to handle init of rC
 */
/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         movddup -128(pB), rC00
         movaps -128(pA), rA0
         movaps rC00, rC10

         mulps rA0, rC00
         movaps rC10, rC20
         movaps -112(pA), rA1

         mulps rA1, rC10
         movddup -120(pB), rC01
         movaps -96(pA), rA2

         mulps rA2, rC20
         movapd rC01, rC11
         movddup -112(pB), rC02

         mulps rA0, rC01
         movaps rC11, rC21
         movddup -104(pB), rB0

         movaps rC02, rC12
         mulps rA1, rC11
         movddup -96(pB), rB1

         mulps rA2, rC21
         movaps rC12, rC22
         prefC(-128(pC))

         mulps rA0, rC02
         movaps -80(pA), rA0
         prefC((pC))

         mulps rA1, rC12
         movaps -64(pA), rA1
         prefB(-128(pfB))

         mulps rA2, rC22
         prefB((pfB))
         add incPF, pfB

         #if KB > 1
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps -48(pA), rA2
               prefetcht0 -88+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup -88(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 2
               movddup -80(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 2
               movaps -32(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 2
               movaps -16(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 2
               movddup -72(pB), rB1
            #endif

         #endif
         #if KB > 2
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 0(pA), rA2
               prefetcht0 -64(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup -64(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 3
               movddup -56(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 3
               movaps 16(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 3
               movaps 32(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 3
               movddup -48(pB), rB1
            #endif

         #endif
         #if KB > 3
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 48(pA), rA2
               prefetcht0 -40+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup -40(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 4
               movddup -32(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 4
               movaps 64(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 4
               movaps 80(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 4
               movddup -24(pB), rB1
            #endif

         #endif
         #if KB > 4
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 96(pA), rA2
               prefetcht0 -16(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup -16(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 5
               movddup -8(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 5
               movaps 112(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 5
               movaps 128(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 5
               movddup 0(pB), rB1
            #endif

         #endif
         #if KB > 5
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 144(pA), rA2
               prefetcht0 8+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 8(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 6
               movddup 16(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 6
               movaps 160(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 6
               movaps 176(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 6
               movddup 24(pB), rB1
            #endif

         #endif
         #if KB > 6
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 192(pA), rA2
               prefetcht0 32(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 32(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 7
               movddup 40(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 7
               movaps 208(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 7
               movaps 224(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 7
               movddup 48(pB), rB1
            #endif

         #endif
         #if KB > 7
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 240(pA), rA2
               prefetcht0 56+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 56(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 8
               movddup 64(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 8
               movaps 256(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 8
               movaps 272(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 8
               movddup 72(pB), rB1
            #endif

         #endif
         #if KB > 8
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 288(pA), rA2
               prefetcht0 80(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 80(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 9
               movddup 88(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 9
               movaps 304(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 9
               movaps 320(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 9
               movddup 96(pB), rB1
            #endif

         #endif
         #if KB > 9
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 336(pA), rA2
               prefetcht0 104+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 104(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 10
               movddup 112(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 10
               movaps 352(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 10
               movaps 368(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 10
               movddup 120(pB), rB1
            #endif

         #endif
         #if KB > 10
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 384(pA), rA2
               prefetcht0 128(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 128(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 11
               movddup 136(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 11
               movaps 400(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 11
               movaps 416(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 11
               movddup 144(pB), rB1
            #endif

         #endif
         #if KB > 11
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 432(pA), rA2
               prefetcht0 152+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 152(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 12
               movddup 160(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 12
               movaps 448(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 12
               movaps 464(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 12
               movddup 168(pB), rB1
            #endif

         #endif
         #if KB > 12
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 480(pA), rA2
               prefetcht0 176(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 176(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 13
               movddup 184(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 13
               movaps 496(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 13
               movaps 512(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 13
               movddup 192(pB), rB1
            #endif

         #endif
         #if KB > 13
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 528(pA), rA2
               prefetcht0 200+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 200(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 14
               movddup 208(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 14
               movaps 544(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 14
               movaps 560(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 14
               movddup 216(pB), rB1
            #endif

         #endif
         #if KB > 14
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 576(pA), rA2
               prefetcht0 224(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 224(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 15
               movddup 232(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 15
               movaps 592(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 15
               movaps 608(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 15
               movddup 240(pB), rB1
            #endif

         #endif
         #if KB > 15
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 624(pA), rA2
               prefetcht0 248+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 248(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 16
               movddup 256(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 16
               movaps 640(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 16
               movaps 656(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 16
               movddup 264(pB), rB1
            #endif

         #endif
         #if KB > 16
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 672(pA), rA2
               prefetcht0 272(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 272(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 17
               movddup 280(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 17
               movaps 688(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 17
               movaps 704(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 17
               movddup 288(pB), rB1
            #endif

         #endif
         #if KB > 17
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 720(pA), rA2
               prefetcht0 296+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 296(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 18
               movddup 304(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 18
               movaps 736(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 18
               movaps 752(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 18
               movddup 312(pB), rB1
            #endif

         #endif
         #if KB > 18
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 768(pA), rA2
               prefetcht0 320(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 320(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 19
               movddup 328(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 19
               movaps 784(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 19
               movaps 800(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 19
               movddup 336(pB), rB1
            #endif

         #endif
         #if KB > 19
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 816(pA), rA2
               prefetcht0 344+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 344(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 20
               movddup 352(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 20
               movaps 832(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 20
               movaps 848(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 20
               movddup 360(pB), rB1
            #endif

         #endif
         #if KB > 20
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 864(pA), rA2
               prefetcht0 368(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 368(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 21
               movddup 376(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 21
               movaps 880(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 21
               movaps 896(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 21
               movddup 384(pB), rB1
            #endif

         #endif
         #if KB > 21
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 912(pA), rA2
               prefetcht0 392+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 392(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 22
               movddup 400(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 22
               movaps 928(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 22
               movaps 944(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 22
               movddup 408(pB), rB1
            #endif

         #endif
         #if KB > 22
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 960(pA), rA2
               prefetcht0 416(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 416(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 23
               movddup 424(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 23
               movaps 976(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 23
               movaps 992(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 23
               movddup 432(pB), rB1
            #endif

         #endif
         #if KB > 23
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1008(pA), rA2
               prefetcht0 440+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 440(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 24
               movddup 448(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 24
               movaps 1024(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 24
               movaps 1040(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 24
               movddup 456(pB), rB1
            #endif

         #endif
         #if KB > 24
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1056(pA), rA2
               prefetcht0 464(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 464(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 25
               movddup 472(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 25
               movaps 1072(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 25
               movaps 1088(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 25
               movddup 480(pB), rB1
            #endif

         #endif
         #if KB > 25
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1104(pA), rA2
               prefetcht0 488+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 488(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 26
               movddup 496(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 26
               movaps 1120(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 26
               movaps 1136(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 26
               movddup 504(pB), rB1
            #endif

         #endif
         #if KB > 26
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1152(pA), rA2
               prefetcht0 512(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 512(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 27
               movddup 520(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 27
               movaps 1168(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 27
               movaps 1184(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 27
               movddup 528(pB), rB1
            #endif

         #endif
         #if KB > 27
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1200(pA), rA2
               prefetcht0 536+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 536(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 28
               movddup 544(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 28
               movaps 1216(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 28
               movaps 1232(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 28
               movddup 552(pB), rB1
            #endif

         #endif
         #if KB > 28
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1248(pA), rA2
               prefetcht0 560(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 560(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 29
               movddup 568(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 29
               movaps 1264(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 29
               movaps 1280(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 29
               movddup 576(pB), rB1
            #endif

         #endif
         #if KB > 29
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1296(pA), rA2
               prefetcht0 584+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 584(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 30
               movddup 592(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 30
               movaps 1312(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 30
               movaps 1328(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 30
               movddup 600(pB), rB1
            #endif

         #endif
         #if KB > 30
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1344(pA), rA2
               prefetcht0 608(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 608(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 31
               movddup 616(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 31
               movaps 1360(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 31
               movaps 1376(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 31
               movddup 624(pB), rB1
            #endif

         #endif
         #if KB > 31
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1392(pA), rA2
               prefetcht0 632+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 632(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 32
               movddup 640(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 32
               movaps 1408(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 32
               movaps 1424(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 32
               movddup 648(pB), rB1
            #endif

         #endif
         #if KB > 32
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1440(pA), rA2
               prefetcht0 656(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 656(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 33
               movddup 664(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 33
               movaps 1456(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 33
               movaps 1472(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 33
               movddup 672(pB), rB1
            #endif

         #endif
         #if KB > 33
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1488(pA), rA2
               prefetcht0 680+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 680(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 34
               movddup 688(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 34
               movaps 1504(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 34
               movaps 1520(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 34
               movddup 696(pB), rB1
            #endif

         #endif
         #if KB > 34
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1536(pA), rA2
               prefetcht0 704(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 704(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 35
               movddup 712(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 35
               movaps 1552(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 35
               movaps 1568(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 35
               movddup 720(pB), rB1
            #endif

         #endif
         #if KB > 35
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1584(pA), rA2
               prefetcht0 728+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 728(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 36
               movddup 736(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 36
               movaps 1600(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 36
               movaps 1616(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 36
               movddup 744(pB), rB1
            #endif

         #endif
         #if KB > 36
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1632(pA), rA2
               prefetcht0 752(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 752(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 37
               movddup 760(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 37
               movaps 1648(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 37
               movaps 1664(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 37
               movddup 768(pB), rB1
            #endif

         #endif
         #if KB > 37
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1680(pA), rA2
               prefetcht0 776+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 776(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 38
               movddup 784(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 38
               movaps 1696(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 38
               movaps 1712(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 38
               movddup 792(pB), rB1
            #endif

         #endif
         #if KB > 38
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1728(pA), rA2
               prefetcht0 800(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 800(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 39
               movddup 808(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 39
               movaps 1744(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 39
               movaps 1760(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 39
               movddup 816(pB), rB1
            #endif

         #endif
         #if KB > 39
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1776(pA), rA2
               prefetcht0 824+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 824(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 40
               movddup 832(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 40
               movaps 1792(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 40
               movaps 1808(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 40
               movddup 840(pB), rB1
            #endif

         #endif
         #if KB > 40
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1824(pA), rA2
               prefetcht0 848(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 848(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 41
               movddup 856(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 41
               movaps 1840(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 41
               movaps 1856(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 41
               movddup 864(pB), rB1
            #endif

         #endif
         #if KB > 41
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1872(pA), rA2
               prefetcht0 872+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 872(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 42
               movddup 880(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 42
               movaps 1888(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 42
               movaps 1904(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 42
               movddup 888(pB), rB1
            #endif

         #endif
         #if KB > 42
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1920(pA), rA2
               prefetcht0 896(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 896(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 43
               movddup 904(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 43
               movaps 1936(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 43
               movaps 1952(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 43
               movddup 912(pB), rB1
            #endif

         #endif
         #if KB > 43
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1968(pA), rA2
               prefetcht0 920+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 920(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 44
               movddup 928(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 44
               movaps 1984(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 44
               movaps 2000(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 44
               movddup 936(pB), rB1
            #endif

         #endif
         #if KB > 44
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2016(pA), rA2
               prefetcht0 944(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 944(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 45
               movddup 952(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 45
               movaps 2032(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 45
               movaps 2048(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 45
               movddup 960(pB), rB1
            #endif

         #endif
         #if KB > 45
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2064(pA), rA2
               prefetcht0 968+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 968(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 46
               movddup 976(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 46
               movaps 2080(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 46
               movaps 2096(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 46
               movddup 984(pB), rB1
            #endif

         #endif
         #if KB > 46
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2112(pA), rA2
               prefetcht0 992(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 992(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 47
               movddup 1000(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 47
               movaps 2128(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 47
               movaps 2144(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 47
               movddup 1008(pB), rB1
            #endif

         #endif
         #if KB > 47
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2160(pA), rA2
               prefetcht0 1016+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1016(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 48
               movddup 1024(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 48
               movaps 2176(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 48
               movaps 2192(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 48
               movddup 1032(pB), rB1
            #endif

         #endif
         #if KB > 48
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2208(pA), rA2
               prefetcht0 1040(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1040(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 49
               movddup 1048(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 49
               movaps 2224(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 49
               movaps 2240(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 49
               movddup 1056(pB), rB1
            #endif

         #endif
         #if KB > 49
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2256(pA), rA2
               prefetcht0 1064+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1064(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 50
               movddup 1072(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 50
               movaps 2272(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 50
               movaps 2288(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 50
               movddup 1080(pB), rB1
            #endif

         #endif
         #if KB > 50
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2304(pA), rA2
               prefetcht0 1088(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1088(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 51
               movddup 1096(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 51
               movaps 2320(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 51
               movaps 2336(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 51
               movddup 1104(pB), rB1
            #endif

         #endif
         #if KB > 51
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2352(pA), rA2
               prefetcht0 1112+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1112(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 52
               movddup 1120(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 52
               movaps 2368(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 52
               movaps 2384(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 52
               movddup 1128(pB), rB1
            #endif

         #endif
         #if KB > 52
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2400(pA), rA2
               prefetcht0 1136(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1136(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 53
               movddup 1144(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 53
               movaps 2416(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 53
               movaps 2432(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 53
               movddup 1152(pB), rB1
            #endif

         #endif
         #if KB > 53
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2448(pA), rA2
               prefetcht0 1160+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1160(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 54
               movddup 1168(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 54
               movaps 2464(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 54
               movaps 2480(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 54
               movddup 1176(pB), rB1
            #endif

         #endif
         #if KB > 54
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2496(pA), rA2
               prefetcht0 1184(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1184(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 55
               movddup 1192(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 55
               movaps 2512(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 55
               movaps 2528(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 55
               movddup 1200(pB), rB1
            #endif

         #endif
         #if KB > 55
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2544(pA), rA2
               prefetcht0 1208+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1208(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 56
               movddup 1216(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 56
               movaps 2560(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 56
               movaps 2576(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 56
               movddup 1224(pB), rB1
            #endif

         #endif
         #if KB > 56
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2592(pA), rA2
               prefetcht0 1232(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1232(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 57
               movddup 1240(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 57
               movaps 2608(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 57
               movaps 2624(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 57
               movddup 1248(pB), rB1
            #endif

         #endif
         #if KB > 57
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2640(pA), rA2
               prefetcht0 1256+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1256(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 58
               movddup 1264(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 58
               movaps 2656(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 58
               movaps 2672(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 58
               movddup 1272(pB), rB1
            #endif

         #endif
         #if KB > 58
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2688(pA), rA2
               prefetcht0 1280(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1280(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 59
               movddup 1288(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 59
               movaps 2704(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 59
               movaps 2720(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 59
               movddup 1296(pB), rB1
            #endif

         #endif
         #if KB > 59
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2736(pA), rA2
               prefetcht0 1304+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1304(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 60
               movddup 1312(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 60
               movaps 2752(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 60
               movaps 2768(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 60
               movddup 1320(pB), rB1
            #endif

         #endif
         #if KB > 60
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2784(pA), rA2
               prefetcht0 1328(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1328(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 61
               movddup 1336(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 61
               movaps 2800(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 61
               movaps 2816(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 61
               movddup 1344(pB), rB1
            #endif

         #endif
         #if KB > 61
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2832(pA), rA2
               prefetcht0 1352+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1352(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 62
               movddup 1360(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 62
               movaps 2848(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 62
               movaps 2864(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 62
               movddup 1368(pB), rB1
            #endif

         #endif
         #if KB > 62
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2880(pA), rA2
               prefetcht0 1376(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1376(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 63
               movddup 1384(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 63
               movaps 2896(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 63
               movaps 2912(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 63
               movddup 1392(pB), rB1
            #endif

         #endif
         #if KB > 63
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2928(pA), rA2
               prefetcht0 1400+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1400(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 64
               movddup 1408(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 64
               movaps 2944(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 64
               movaps 2960(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 64
               movddup 1416(pB), rB1
            #endif

         #endif
         #if KB > 64
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2976(pA), rA2
               prefetcht0 1424(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1424(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 65
               movddup 1432(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 65
               movaps 2992(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 65
               movaps 3008(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 65
               movddup 1440(pB), rB1
            #endif

         #endif
         #if KB > 65
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3024(pA), rA2
               prefetcht0 1448+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1448(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 66
               movddup 1456(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 66
               movaps 3040(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 66
               movaps 3056(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 66
               movddup 1464(pB), rB1
            #endif

         #endif
         #if KB > 66
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3072(pA), rA2
               prefetcht0 1472(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1472(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 67
               movddup 1480(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 67
               movaps 3088(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 67
               movaps 3104(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 67
               movddup 1488(pB), rB1
            #endif

         #endif
         #if KB > 67
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3120(pA), rA2
               prefetcht0 1496+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1496(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 68
               movddup 1504(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 68
               movaps 3136(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 68
               movaps 3152(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 68
               movddup 1512(pB), rB1
            #endif

         #endif
         #if KB > 68
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3168(pA), rA2
               prefetcht0 1520(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1520(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 69
               movddup 1528(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 69
               movaps 3184(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 69
               movaps 3200(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 69
               movddup 1536(pB), rB1
            #endif

         #endif
         #if KB > 69
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3216(pA), rA2
               prefetcht0 1544+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1544(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 70
               movddup 1552(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 70
               movaps 3232(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 70
               movaps 3248(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 70
               movddup 1560(pB), rB1
            #endif

         #endif
         #if KB > 70
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3264(pA), rA2
               prefetcht0 1568(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1568(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 71
               movddup 1576(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 71
               movaps 3280(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 71
               movaps 3296(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 71
               movddup 1584(pB), rB1
            #endif

         #endif
         #if KB > 71
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3312(pA), rA2
               prefetcht0 1592+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1592(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 72
               movddup 1600(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 72
               movaps 3328(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 72
               movaps 3344(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 72
               movddup 1608(pB), rB1
            #endif

         #endif
         #if KB > 72
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3360(pA), rA2
               prefetcht0 1616(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1616(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 73
               movddup 1624(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 73
               movaps 3376(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 73
               movaps 3392(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 73
               movddup 1632(pB), rB1
            #endif

         #endif
         #if KB > 73
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3408(pA), rA2
               prefetcht0 1640+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1640(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 74
               movddup 1648(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 74
               movaps 3424(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 74
               movaps 3440(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 74
               movddup 1656(pB), rB1
            #endif

         #endif
         #if KB > 74
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3456(pA), rA2
               prefetcht0 1664(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1664(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 75
               movddup 1672(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 75
               movaps 3472(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 75
               movaps 3488(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 75
               movddup 1680(pB), rB1
            #endif

         #endif
         #if KB > 75
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3504(pA), rA2
               prefetcht0 1688+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1688(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 76
               movddup 1696(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 76
               movaps 3520(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 76
               movaps 3536(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 76
               movddup 1704(pB), rB1
            #endif

         #endif
         #if KB > 76
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3552(pA), rA2
               prefetcht0 1712(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1712(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 77
               movddup 1720(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 77
               movaps 3568(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 77
               movaps 3584(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 77
               movddup 1728(pB), rB1
            #endif

         #endif
         #if KB > 77
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3600(pA), rA2
               prefetcht0 1736+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1736(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 78
               movddup 1744(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 78
               movaps 3616(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 78
               movaps 3632(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 78
               movddup 1752(pB), rB1
            #endif

         #endif
         #if KB > 78
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3648(pA), rA2
               prefetcht0 1760(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1760(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 79
               movddup 1768(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 79
               movaps 3664(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 79
               movaps 3680(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 79
               movddup 1776(pB), rB1
            #endif

         #endif
         #if KB > 79
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3696(pA), rA2
               prefetcht0 1784+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1784(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 80
               movddup 1792(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 80
               movaps 3712(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 80
               movaps 3728(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 80
               movddup 1800(pB), rB1
            #endif

         #endif
         #if KB > 80
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3744(pA), rA2
               prefetcht0 1808(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1808(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 81
               movddup 1816(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 81
               movaps 3760(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 81
               movaps 3776(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 81
               movddup 1824(pB), rB1
            #endif

         #endif
         #if KB > 81
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3792(pA), rA2
               prefetcht0 1832+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1832(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 82
               movddup 1840(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 82
               movaps 3808(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 82
               movaps 3824(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 82
               movddup 1848(pB), rB1
            #endif

         #endif
         #if KB > 82
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3840(pA), rA2
               prefetcht0 1856(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1856(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 83
               movddup 1864(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 83
               movaps 3856(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 83
               movaps 3872(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 83
               movddup 1872(pB), rB1
            #endif

         #endif
         #if KB > 83
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3888(pA), rA2
               prefetcht0 1880+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1880(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 84
               movddup 1888(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 84
               movaps 3904(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 84
               movaps 3920(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 84
               movddup 1896(pB), rB1
            #endif

         #endif
         #if KB > 84
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3936(pA), rA2
               prefetcht0 1904(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1904(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 85
               movddup 1912(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 85
               movaps 3952(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 85
               movaps 3968(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 85
               movddup 1920(pB), rB1
            #endif

         #endif
         #if KB > 85
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3984(pA), rA2
               prefetcht0 1928+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1928(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 86
               movddup 1936(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 86
               movaps 4000(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 86
               movaps 4016(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 86
               movddup 1944(pB), rB1
            #endif

         #endif
         #if KB > 86
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4032(pA), rA2
               prefetcht0 1952(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1952(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 87
               movddup 1960(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 87
               movaps 4048(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 87
               movaps 4064(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 87
               movddup 1968(pB), rB1
            #endif

         #endif
         #if KB > 87
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4080(pA), rA2
               prefetcht0 1976+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 1976(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 88
               movddup 1984(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 88
               movaps 4096(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 88
               movaps 4112(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 88
               movddup 1992(pB), rB1
            #endif

         #endif
         #if KB > 88
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4128(pA), rA2
               prefetcht0 2000(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2000(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 89
               movddup 2008(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 89
               movaps 4144(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 89
               movaps 4160(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 89
               movddup 2016(pB), rB1
            #endif

         #endif
         #if KB > 89
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4176(pA), rA2
               prefetcht0 2024+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2024(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 90
               movddup 2032(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 90
               movaps 4192(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 90
               movaps 4208(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 90
               movddup 2040(pB), rB1
            #endif

         #endif
         #if KB > 90
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4224(pA), rA2
               prefetcht0 2048(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2048(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 91
               movddup 2056(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 91
               movaps 4240(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 91
               movaps 4256(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 91
               movddup 2064(pB), rB1
            #endif

         #endif
         #if KB > 91
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4272(pA), rA2
               prefetcht0 2072+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2072(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 92
               movddup 2080(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 92
               movaps 4288(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 92
               movaps 4304(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 92
               movddup 2088(pB), rB1
            #endif

         #endif
         #if KB > 92
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4320(pA), rA2
               prefetcht0 2096(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2096(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 93
               movddup 2104(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 93
               movaps 4336(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 93
               movaps 4352(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 93
               movddup 2112(pB), rB1
            #endif

         #endif
         #if KB > 93
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4368(pA), rA2
               prefetcht0 2120+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2120(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 94
               movddup 2128(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 94
               movaps 4384(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 94
               movaps 4400(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 94
               movddup 2136(pB), rB1
            #endif

         #endif
         #if KB > 94
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4416(pA), rA2
               prefetcht0 2144(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2144(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 95
               movddup 2152(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 95
               movaps 4432(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 95
               movaps 4448(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 95
               movddup 2160(pB), rB1
            #endif

         #endif
         #if KB > 95
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4464(pA), rA2
               prefetcht0 2168+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2168(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 96
               movddup 2176(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 96
               movaps 4480(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 96
               movaps 4496(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 96
               movddup 2184(pB), rB1
            #endif

         #endif
         #if KB > 96
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4512(pA), rA2
               prefetcht0 2192(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2192(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 97
               movddup 2200(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 97
               movaps 4528(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 97
               movaps 4544(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 97
               movddup 2208(pB), rB1
            #endif

         #endif
         #if KB > 97
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4560(pA), rA2
               prefetcht0 2216+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2216(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 98
               movddup 2224(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 98
               movaps 4576(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 98
               movaps 4592(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 98
               movddup 2232(pB), rB1
            #endif

         #endif
         #if KB > 98
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4608(pA), rA2
               prefetcht0 2240(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2240(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 99
               movddup 2248(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 99
               movaps 4624(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 99
               movaps 4640(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 99
               movddup 2256(pB), rB1
            #endif

         #endif
         #if KB > 99
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4656(pA), rA2
               prefetcht0 2264+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2264(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 100
               movddup 2272(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 100
               movaps 4672(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 100
               movaps 4688(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 100
               movddup 2280(pB), rB1
            #endif

         #endif
         #if KB > 100
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4704(pA), rA2
               prefetcht0 2288(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2288(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 101
               movddup 2296(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 101
               movaps 4720(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 101
               movaps 4736(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 101
               movddup 2304(pB), rB1
            #endif

         #endif
         #if KB > 101
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4752(pA), rA2
               prefetcht0 2312+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2312(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 102
               movddup 2320(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 102
               movaps 4768(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 102
               movaps 4784(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 102
               movddup 2328(pB), rB1
            #endif

         #endif
         #if KB > 102
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4800(pA), rA2
               prefetcht0 2336(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2336(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 103
               movddup 2344(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 103
               movaps 4816(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 103
               movaps 4832(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 103
               movddup 2352(pB), rB1
            #endif

         #endif
         #if KB > 103
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4848(pA), rA2
               prefetcht0 2360+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2360(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 104
               movddup 2368(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 104
               movaps 4864(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 104
               movaps 4880(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 104
               movddup 2376(pB), rB1
            #endif

         #endif
         #if KB > 104
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4896(pA), rA2
               prefetcht0 2384(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2384(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 105
               movddup 2392(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 105
               movaps 4912(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 105
               movaps 4928(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 105
               movddup 2400(pB), rB1
            #endif

         #endif
         #if KB > 105
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4944(pA), rA2
               prefetcht0 2408+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2408(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 106
               movddup 2416(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 106
               movaps 4960(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 106
               movaps 4976(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 106
               movddup 2424(pB), rB1
            #endif

         #endif
         #if KB > 106
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4992(pA), rA2
               prefetcht0 2432(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2432(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 107
               movddup 2440(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 107
               movaps 5008(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 107
               movaps 5024(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 107
               movddup 2448(pB), rB1
            #endif

         #endif
         #if KB > 107
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5040(pA), rA2
               prefetcht0 2456+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2456(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 108
               movddup 2464(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 108
               movaps 5056(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 108
               movaps 5072(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 108
               movddup 2472(pB), rB1
            #endif

         #endif
         #if KB > 108
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5088(pA), rA2
               prefetcht0 2480(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2480(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 109
               movddup 2488(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 109
               movaps 5104(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 109
               movaps 5120(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 109
               movddup 2496(pB), rB1
            #endif

         #endif
         #if KB > 109
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5136(pA), rA2
               prefetcht0 2504+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2504(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 110
               movddup 2512(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 110
               movaps 5152(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 110
               movaps 5168(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 110
               movddup 2520(pB), rB1
            #endif

         #endif
         #if KB > 110
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5184(pA), rA2
               prefetcht0 2528(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2528(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 111
               movddup 2536(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 111
               movaps 5200(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 111
               movaps 5216(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 111
               movddup 2544(pB), rB1
            #endif

         #endif
         #if KB > 111
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5232(pA), rA2
               prefetcht0 2552+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2552(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 112
               movddup 2560(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 112
               movaps 5248(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 112
               movaps 5264(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 112
               movddup 2568(pB), rB1
            #endif

         #endif
         #if KB > 112
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5280(pA), rA2
               prefetcht0 2576(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2576(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 113
               movddup 2584(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 113
               movaps 5296(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 113
               movaps 5312(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 113
               movddup 2592(pB), rB1
            #endif

         #endif
         #if KB > 113
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5328(pA), rA2
               prefetcht0 2600+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2600(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 114
               movddup 2608(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 114
               movaps 5344(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 114
               movaps 5360(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 114
               movddup 2616(pB), rB1
            #endif

         #endif
         #if KB > 114
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5376(pA), rA2
               prefetcht0 2624(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2624(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 115
               movddup 2632(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 115
               movaps 5392(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 115
               movaps 5408(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 115
               movddup 2640(pB), rB1
            #endif

         #endif
         #if KB > 115
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5424(pA), rA2
               prefetcht0 2648+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2648(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 116
               movddup 2656(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 116
               movaps 5440(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 116
               movaps 5456(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 116
               movddup 2664(pB), rB1
            #endif

         #endif
         #if KB > 116
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5472(pA), rA2
               prefetcht0 2672(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2672(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 117
               movddup 2680(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 117
               movaps 5488(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 117
               movaps 5504(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 117
               movddup 2688(pB), rB1
            #endif

         #endif
         #if KB > 117
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5520(pA), rA2
               prefetcht0 2696+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2696(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 118
               movddup 2704(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 118
               movaps 5536(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 118
               movaps 5552(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 118
               movddup 2712(pB), rB1
            #endif

         #endif
         #if KB > 118
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5568(pA), rA2
               prefetcht0 2720(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2720(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 119
               movddup 2728(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 119
               movaps 5584(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 119
               movaps 5600(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 119
               movddup 2736(pB), rB1
            #endif

         #endif
         #if KB > 119
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5616(pA), rA2
               prefetcht0 2744+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2744(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 120
               movddup 2752(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 120
               movaps 5632(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 120
               movaps 5648(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 120
               movddup 2760(pB), rB1
            #endif

         #endif
         #if KB > 120
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5664(pA), rA2
               prefetcht0 2768(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2768(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 121
               movddup 2776(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 121
               movaps 5680(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 121
               movaps 5696(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 121
               movddup 2784(pB), rB1
            #endif

         #endif
         #if KB > 121
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5712(pA), rA2
               prefetcht0 2792+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2792(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 122
               movddup 2800(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 122
               movaps 5728(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 122
               movaps 5744(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 122
               movddup 2808(pB), rB1
            #endif

         #endif
         #if KB > 122
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5760(pA), rA2
               prefetcht0 2816(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2816(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 123
               movddup 2824(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 123
               movaps 5776(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 123
               movaps 5792(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 123
               movddup 2832(pB), rB1
            #endif

         #endif
         #if KB > 123
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5808(pA), rA2
               prefetcht0 2840+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2840(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 124
               movddup 2848(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 124
               movaps 5824(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 124
               movaps 5840(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 124
               movddup 2856(pB), rB1
            #endif

         #endif
         #if KB > 124
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5856(pA), rA2
               prefetcht0 2864(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2864(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 125
               movddup 2872(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 125
               movaps 5872(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 125
               movaps 5888(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 125
               movddup 2880(pB), rB1
            #endif

         #endif
         #if KB > 125
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5904(pA), rA2
               prefetcht0 2888+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2888(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 126
               movddup 2896(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 126
               movaps 5920(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 126
               movaps 5936(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 126
               movddup 2904(pB), rB1
            #endif

         #endif
         #if KB > 126
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5952(pA), rA2
               prefetcht0 2912(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2912(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 127
               movddup 2920(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 127
               movaps 5968(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 127
               movaps 5984(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 127
               movddup 2928(pB), rB1
            #endif

         #endif
         #if KB > 127
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6000(pA), rA2
               prefetcht0 2936+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2936(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 128
               movddup 2944(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 128
               movaps 6016(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 128
               movaps 6032(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 128
               movddup 2952(pB), rB1
            #endif

         #endif
         #if KB > 128
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6048(pA), rA2
               prefetcht0 2960(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2960(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 129
               movddup 2968(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 129
               movaps 6064(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 129
               movaps 6080(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 129
               movddup 2976(pB), rB1
            #endif

         #endif
         #if KB > 129
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6096(pA), rA2
               prefetcht0 2984+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 2984(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 130
               movddup 2992(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 130
               movaps 6112(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 130
               movaps 6128(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 130
               movddup 3000(pB), rB1
            #endif

         #endif
         #if KB > 130
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6144(pA), rA2
               prefetcht0 3008(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3008(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 131
               movddup 3016(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 131
               movaps 6160(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 131
               movaps 6176(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 131
               movddup 3024(pB), rB1
            #endif

         #endif
         #if KB > 131
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6192(pA), rA2
               prefetcht0 3032+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3032(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 132
               movddup 3040(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 132
               movaps 6208(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 132
               movaps 6224(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 132
               movddup 3048(pB), rB1
            #endif

         #endif
         #if KB > 132
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6240(pA), rA2
               prefetcht0 3056(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3056(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 133
               movddup 3064(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 133
               movaps 6256(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 133
               movaps 6272(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 133
               movddup 3072(pB), rB1
            #endif

         #endif
         #if KB > 133
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6288(pA), rA2
               prefetcht0 3080+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3080(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 134
               movddup 3088(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 134
               movaps 6304(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 134
               movaps 6320(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 134
               movddup 3096(pB), rB1
            #endif

         #endif
         #if KB > 134
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6336(pA), rA2
               prefetcht0 3104(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3104(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 135
               movddup 3112(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 135
               movaps 6352(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 135
               movaps 6368(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 135
               movddup 3120(pB), rB1
            #endif

         #endif
         #if KB > 135
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6384(pA), rA2
               prefetcht0 3128+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3128(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 136
               movddup 3136(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 136
               movaps 6400(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 136
               movaps 6416(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 136
               movddup 3144(pB), rB1
            #endif

         #endif
         #if KB > 136
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6432(pA), rA2
               prefetcht0 3152(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3152(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 137
               movddup 3160(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 137
               movaps 6448(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 137
               movaps 6464(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 137
               movddup 3168(pB), rB1
            #endif

         #endif
         #if KB > 137
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6480(pA), rA2
               prefetcht0 3176+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3176(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 138
               movddup 3184(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 138
               movaps 6496(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 138
               movaps 6512(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 138
               movddup 3192(pB), rB1
            #endif

         #endif
         #if KB > 138
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6528(pA), rA2
               prefetcht0 3200(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3200(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 139
               movddup 3208(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 139
               movaps 6544(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 139
               movaps 6560(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 139
               movddup 3216(pB), rB1
            #endif

         #endif
         #if KB > 139
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6576(pA), rA2
               prefetcht0 3224+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3224(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 140
               movddup 3232(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 140
               movaps 6592(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 140
               movaps 6608(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 140
               movddup 3240(pB), rB1
            #endif

         #endif
         #if KB > 140
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6624(pA), rA2
               prefetcht0 3248(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3248(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 141
               movddup 3256(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 141
               movaps 6640(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 141
               movaps 6656(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 141
               movddup 3264(pB), rB1
            #endif

         #endif
         #if KB > 141
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6672(pA), rA2
               prefetcht0 3272+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3272(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 142
               movddup 3280(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 142
               movaps 6688(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 142
               movaps 6704(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 142
               movddup 3288(pB), rB1
            #endif

         #endif
         #if KB > 142
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6720(pA), rA2
               prefetcht0 3296(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3296(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 143
               movddup 3304(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 143
               movaps 6736(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 143
               movaps 6752(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 143
               movddup 3312(pB), rB1
            #endif

         #endif
         #if KB > 143
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6768(pA), rA2
               prefetcht0 3320+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3320(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 144
               movddup 3328(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 144
               movaps 6784(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 144
               movaps 6800(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 144
               movddup 3336(pB), rB1
            #endif

         #endif
         #if KB > 144
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6816(pA), rA2
               prefetcht0 3344(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3344(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 145
               movddup 3352(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 145
               movaps 6832(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 145
               movaps 6848(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 145
               movddup 3360(pB), rB1
            #endif

         #endif
         #if KB > 145
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6864(pA), rA2
               prefetcht0 3368+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3368(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 146
               movddup 3376(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 146
               movaps 6880(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 146
               movaps 6896(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 146
               movddup 3384(pB), rB1
            #endif

         #endif
         #if KB > 146
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6912(pA), rA2
               prefetcht0 3392(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3392(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 147
               movddup 3400(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 147
               movaps 6928(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 147
               movaps 6944(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 147
               movddup 3408(pB), rB1
            #endif

         #endif
         #if KB > 147
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6960(pA), rA2
               prefetcht0 3416+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3416(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 148
               movddup 3424(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 148
               movaps 6976(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 148
               movaps 6992(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 148
               movddup 3432(pB), rB1
            #endif

         #endif
         #if KB > 148
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7008(pA), rA2
               prefetcht0 3440(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3440(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 149
               movddup 3448(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 149
               movaps 7024(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 149
               movaps 7040(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 149
               movddup 3456(pB), rB1
            #endif

         #endif
         #if KB > 149
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7056(pA), rA2
               prefetcht0 3464+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3464(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 150
               movddup 3472(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 150
               movaps 7072(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 150
               movaps 7088(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 150
               movddup 3480(pB), rB1
            #endif

         #endif
         #if KB > 150
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7104(pA), rA2
               prefetcht0 3488(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3488(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 151
               movddup 3496(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 151
               movaps 7120(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 151
               movaps 7136(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 151
               movddup 3504(pB), rB1
            #endif

         #endif
         #if KB > 151
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7152(pA), rA2
               prefetcht0 3512+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3512(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 152
               movddup 3520(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 152
               movaps 7168(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 152
               movaps 7184(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 152
               movddup 3528(pB), rB1
            #endif

         #endif
         #if KB > 152
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7200(pA), rA2
               prefetcht0 3536(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3536(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 153
               movddup 3544(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 153
               movaps 7216(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 153
               movaps 7232(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 153
               movddup 3552(pB), rB1
            #endif

         #endif
         #if KB > 153
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7248(pA), rA2
               prefetcht0 3560+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3560(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 154
               movddup 3568(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 154
               movaps 7264(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 154
               movaps 7280(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 154
               movddup 3576(pB), rB1
            #endif

         #endif
         #if KB > 154
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7296(pA), rA2
               prefetcht0 3584(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3584(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 155
               movddup 3592(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 155
               movaps 7312(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 155
               movaps 7328(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 155
               movddup 3600(pB), rB1
            #endif

         #endif
         #if KB > 155
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7344(pA), rA2
               prefetcht0 3608+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3608(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 156
               movddup 3616(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 156
               movaps 7360(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 156
               movaps 7376(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 156
               movddup 3624(pB), rB1
            #endif

         #endif
         #if KB > 156
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7392(pA), rA2
               prefetcht0 3632(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3632(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 157
               movddup 3640(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 157
               movaps 7408(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 157
               movaps 7424(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 157
               movddup 3648(pB), rB1
            #endif

         #endif
         #if KB > 157
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7440(pA), rA2
               prefetcht0 3656+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3656(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 158
               movddup 3664(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 158
               movaps 7456(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 158
               movaps 7472(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 158
               movddup 3672(pB), rB1
            #endif

         #endif
         #if KB > 158
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7488(pA), rA2
               prefetcht0 3680(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3680(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 159
               movddup 3688(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 159
               movaps 7504(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 159
               movaps 7520(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 159
               movddup 3696(pB), rB1
            #endif

         #endif
         #if KB > 159
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7536(pA), rA2
               prefetcht0 3704+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3704(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 160
               movddup 3712(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 160
               movaps 7552(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 160
               movaps 7568(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 160
               movddup 3720(pB), rB1
            #endif

         #endif
         #if KB > 160
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7584(pA), rA2
               prefetcht0 3728(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3728(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 161
               movddup 3736(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 161
               movaps 7600(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 161
               movaps 7616(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 161
               movddup 3744(pB), rB1
            #endif

         #endif
         #if KB > 161
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7632(pA), rA2
               prefetcht0 3752+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3752(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 162
               movddup 3760(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 162
               movaps 7648(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 162
               movaps 7664(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 162
               movddup 3768(pB), rB1
            #endif

         #endif
         #if KB > 162
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7680(pA), rA2
               prefetcht0 3776(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3776(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 163
               movddup 3784(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 163
               movaps 7696(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 163
               movaps 7712(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 163
               movddup 3792(pB), rB1
            #endif

         #endif
         #if KB > 163
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7728(pA), rA2
               prefetcht0 3800+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3800(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 164
               movddup 3808(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 164
               movaps 7744(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 164
               movaps 7760(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 164
               movddup 3816(pB), rB1
            #endif

         #endif
         #if KB > 164
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7776(pA), rA2
               prefetcht0 3824(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3824(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 165
               movddup 3832(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 165
               movaps 7792(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 165
               movaps 7808(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 165
               movddup 3840(pB), rB1
            #endif

         #endif
         #if KB > 165
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7824(pA), rA2
               prefetcht0 3848+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3848(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 166
               movddup 3856(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 166
               movaps 7840(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 166
               movaps 7856(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 166
               movddup 3864(pB), rB1
            #endif

         #endif
         #if KB > 166
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7872(pA), rA2
               prefetcht0 3872(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3872(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 167
               movddup 3880(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 167
               movaps 7888(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 167
               movaps 7904(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 167
               movddup 3888(pB), rB1
            #endif

         #endif
         #if KB > 167
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7920(pA), rA2
               prefetcht0 3896+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3896(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 168
               movddup 3904(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 168
               movaps 7936(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 168
               movaps 7952(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 168
               movddup 3912(pB), rB1
            #endif

         #endif
         #if KB > 168
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7968(pA), rA2
               prefetcht0 3920(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3920(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 169
               movddup 3928(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 169
               movaps 7984(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 169
               movaps 8000(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 169
               movddup 3936(pB), rB1
            #endif

         #endif
         #if KB > 169
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8016(pA), rA2
               prefetcht0 3944+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3944(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 170
               movddup 3952(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 170
               movaps 8032(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 170
               movaps 8048(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 170
               movddup 3960(pB), rB1
            #endif

         #endif
         #if KB > 170
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8064(pA), rA2
               prefetcht0 3968(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3968(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 171
               movddup 3976(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 171
               movaps 8080(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 171
               movaps 8096(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 171
               movddup 3984(pB), rB1
            #endif

         #endif
         #if KB > 171
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8112(pA), rA2
               prefetcht0 3992+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 3992(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 172
               movddup 4000(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 172
               movaps 8128(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 172
               movaps 8144(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 172
               movddup 4008(pB), rB1
            #endif

         #endif
         #if KB > 172
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8160(pA), rA2
               prefetcht0 4016(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4016(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 173
               movddup 4024(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 173
               movaps 8176(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 173
               movaps 8192(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 173
               movddup 4032(pB), rB1
            #endif

         #endif
         #if KB > 173
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8208(pA), rA2
               prefetcht0 4040+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4040(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 174
               movddup 4048(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 174
               movaps 8224(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 174
               movaps 8240(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 174
               movddup 4056(pB), rB1
            #endif

         #endif
         #if KB > 174
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8256(pA), rA2
               prefetcht0 4064(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4064(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 175
               movddup 4072(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 175
               movaps 8272(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 175
               movaps 8288(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 175
               movddup 4080(pB), rB1
            #endif

         #endif
         #if KB > 175
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8304(pA), rA2
               prefetcht0 4088+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4088(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 176
               movddup 4096(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 176
               movaps 8320(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 176
               movaps 8336(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 176
               movddup 4104(pB), rB1
            #endif

         #endif
         #if KB > 176
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8352(pA), rA2
               prefetcht0 4112(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4112(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 177
               movddup 4120(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 177
               movaps 8368(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 177
               movaps 8384(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 177
               movddup 4128(pB), rB1
            #endif

         #endif
         #if KB > 177
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8400(pA), rA2
               prefetcht0 4136+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4136(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 178
               movddup 4144(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 178
               movaps 8416(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 178
               movaps 8432(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 178
               movddup 4152(pB), rB1
            #endif

         #endif
         #if KB > 178
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8448(pA), rA2
               prefetcht0 4160(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4160(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 179
               movddup 4168(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 179
               movaps 8464(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 179
               movaps 8480(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 179
               movddup 4176(pB), rB1
            #endif

         #endif
         #if KB > 179
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8496(pA), rA2
               prefetcht0 4184+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4184(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 180
               movddup 4192(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 180
               movaps 8512(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 180
               movaps 8528(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 180
               movddup 4200(pB), rB1
            #endif

         #endif
         #if KB > 180
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8544(pA), rA2
               prefetcht0 4208(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4208(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 181
               movddup 4216(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 181
               movaps 8560(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 181
               movaps 8576(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 181
               movddup 4224(pB), rB1
            #endif

         #endif
         #if KB > 181
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8592(pA), rA2
               prefetcht0 4232+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4232(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 182
               movddup 4240(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 182
               movaps 8608(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 182
               movaps 8624(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 182
               movddup 4248(pB), rB1
            #endif

         #endif
         #if KB > 182
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8640(pA), rA2
               prefetcht0 4256(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4256(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 183
               movddup 4264(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 183
               movaps 8656(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 183
               movaps 8672(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 183
               movddup 4272(pB), rB1
            #endif

         #endif
         #if KB > 183
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8688(pA), rA2
               prefetcht0 4280+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4280(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 184
               movddup 4288(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 184
               movaps 8704(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 184
               movaps 8720(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 184
               movddup 4296(pB), rB1
            #endif

         #endif
         #if KB > 184
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8736(pA), rA2
               prefetcht0 4304(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4304(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 185
               movddup 4312(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 185
               movaps 8752(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 185
               movaps 8768(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 185
               movddup 4320(pB), rB1
            #endif

         #endif
         #if KB > 185
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8784(pA), rA2
               prefetcht0 4328+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4328(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 186
               movddup 4336(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 186
               movaps 8800(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 186
               movaps 8816(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 186
               movddup 4344(pB), rB1
            #endif

         #endif
         #if KB > 186
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8832(pA), rA2
               prefetcht0 4352(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4352(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 187
               movddup 4360(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 187
               movaps 8848(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 187
               movaps 8864(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 187
               movddup 4368(pB), rB1
            #endif

         #endif
         #if KB > 187
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8880(pA), rA2
               prefetcht0 4376+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4376(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 188
               movddup 4384(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 188
               movaps 8896(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 188
               movaps 8912(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 188
               movddup 4392(pB), rB1
            #endif

         #endif
         #if KB > 188
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8928(pA), rA2
               prefetcht0 4400(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4400(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 189
               movddup 4408(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 189
               movaps 8944(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 189
               movaps 8960(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 189
               movddup 4416(pB), rB1
            #endif

         #endif
         #if KB > 189
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8976(pA), rA2
               prefetcht0 4424+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4424(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 190
               movddup 4432(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 190
               movaps 8992(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 190
               movaps 9008(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 190
               movddup 4440(pB), rB1
            #endif

         #endif
         #if KB > 190
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9024(pA), rA2
               prefetcht0 4448(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4448(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 191
               movddup 4456(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 191
               movaps 9040(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 191
               movaps 9056(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 191
               movddup 4464(pB), rB1
            #endif

         #endif
         #if KB > 191
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9072(pA), rA2
               prefetcht0 4472+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4472(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 192
               movddup 4480(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 192
               movaps 9088(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 192
               movaps 9104(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 192
               movddup 4488(pB), rB1
            #endif

         #endif
         #if KB > 192
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9120(pA), rA2
               prefetcht0 4496(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4496(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 193
               movddup 4504(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 193
               movaps 9136(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 193
               movaps 9152(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 193
               movddup 4512(pB), rB1
            #endif

         #endif
         #if KB > 193
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9168(pA), rA2
               prefetcht0 4520+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4520(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 194
               movddup 4528(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 194
               movaps 9184(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 194
               movaps 9200(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 194
               movddup 4536(pB), rB1
            #endif

         #endif
         #if KB > 194
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9216(pA), rA2
               prefetcht0 4544(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4544(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 195
               movddup 4552(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 195
               movaps 9232(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 195
               movaps 9248(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 195
               movddup 4560(pB), rB1
            #endif

         #endif
         #if KB > 195
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9264(pA), rA2
               prefetcht0 4568+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4568(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 196
               movddup 4576(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 196
               movaps 9280(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 196
               movaps 9296(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 196
               movddup 4584(pB), rB1
            #endif

         #endif
         #if KB > 196
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9312(pA), rA2
               prefetcht0 4592(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4592(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 197
               movddup 4600(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 197
               movaps 9328(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 197
               movaps 9344(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 197
               movddup 4608(pB), rB1
            #endif

         #endif
         #if KB > 197
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9360(pA), rA2
               prefetcht0 4616+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4616(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 198
               movddup 4624(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 198
               movaps 9376(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 198
               movaps 9392(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 198
               movddup 4632(pB), rB1
            #endif

         #endif
         #if KB > 198
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9408(pA), rA2
               prefetcht0 4640(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4640(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 199
               movddup 4648(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 199
               movaps 9424(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 199
               movaps 9440(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 199
               movddup 4656(pB), rB1
            #endif

         #endif
         #if KB > 199
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9456(pA), rA2
               prefetcht0 4664+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4664(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 200
               movddup 4672(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 200
               movaps 9472(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 200
               movaps 9488(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 200
               movddup 4680(pB), rB1
            #endif

         #endif
         #if KB > 200
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9504(pA), rA2
               prefetcht0 4688(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4688(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 201
               movddup 4696(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 201
               movaps 9520(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 201
               movaps 9536(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 201
               movddup 4704(pB), rB1
            #endif

         #endif
         #if KB > 201
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9552(pA), rA2
               prefetcht0 4712+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4712(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 202
               movddup 4720(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 202
               movaps 9568(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 202
               movaps 9584(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 202
               movddup 4728(pB), rB1
            #endif

         #endif
         #if KB > 202
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9600(pA), rA2
               prefetcht0 4736(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4736(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 203
               movddup 4744(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 203
               movaps 9616(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 203
               movaps 9632(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 203
               movddup 4752(pB), rB1
            #endif

         #endif
         #if KB > 203
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9648(pA), rA2
               prefetcht0 4760+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4760(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 204
               movddup 4768(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 204
               movaps 9664(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 204
               movaps 9680(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 204
               movddup 4776(pB), rB1
            #endif

         #endif
         #if KB > 204
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9696(pA), rA2
               prefetcht0 4784(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4784(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 205
               movddup 4792(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 205
               movaps 9712(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 205
               movaps 9728(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 205
               movddup 4800(pB), rB1
            #endif

         #endif
         #if KB > 205
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9744(pA), rA2
               prefetcht0 4808+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4808(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 206
               movddup 4816(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 206
               movaps 9760(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 206
               movaps 9776(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 206
               movddup 4824(pB), rB1
            #endif

         #endif
         #if KB > 206
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9792(pA), rA2
               prefetcht0 4832(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4832(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 207
               movddup 4840(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 207
               movaps 9808(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 207
               movaps 9824(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 207
               movddup 4848(pB), rB1
            #endif

         #endif
         #if KB > 207
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9840(pA), rA2
               prefetcht0 4856+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4856(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 208
               movddup 4864(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 208
               movaps 9856(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 208
               movaps 9872(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 208
               movddup 4872(pB), rB1
            #endif

         #endif
         #if KB > 208
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9888(pA), rA2
               prefetcht0 4880(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4880(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 209
               movddup 4888(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 209
               movaps 9904(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 209
               movaps 9920(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 209
               movddup 4896(pB), rB1
            #endif

         #endif
         #if KB > 209
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9936(pA), rA2
               prefetcht0 4904+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4904(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 210
               movddup 4912(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 210
               movaps 9952(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 210
               movaps 9968(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 210
               movddup 4920(pB), rB1
            #endif

         #endif
         #if KB > 210
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9984(pA), rA2
               prefetcht0 4928(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4928(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 211
               movddup 4936(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 211
               movaps 10000(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 211
               movaps 10016(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 211
               movddup 4944(pB), rB1
            #endif

         #endif
         #if KB > 211
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10032(pA), rA2
               prefetcht0 4952+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4952(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 212
               movddup 4960(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 212
               movaps 10048(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 212
               movaps 10064(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 212
               movddup 4968(pB), rB1
            #endif

         #endif
         #if KB > 212
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10080(pA), rA2
               prefetcht0 4976(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 4976(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 213
               movddup 4984(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 213
               movaps 10096(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 213
               movaps 10112(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 213
               movddup 4992(pB), rB1
            #endif

         #endif
         #if KB > 213
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10128(pA), rA2
               prefetcht0 5000+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5000(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 214
               movddup 5008(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 214
               movaps 10144(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 214
               movaps 10160(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 214
               movddup 5016(pB), rB1
            #endif

         #endif
         #if KB > 214
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10176(pA), rA2
               prefetcht0 5024(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5024(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 215
               movddup 5032(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 215
               movaps 10192(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 215
               movaps 10208(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 215
               movddup 5040(pB), rB1
            #endif

         #endif
         #if KB > 215
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10224(pA), rA2
               prefetcht0 5048+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5048(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 216
               movddup 5056(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 216
               movaps 10240(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 216
               movaps 10256(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 216
               movddup 5064(pB), rB1
            #endif

         #endif
         #if KB > 216
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10272(pA), rA2
               prefetcht0 5072(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5072(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 217
               movddup 5080(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 217
               movaps 10288(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 217
               movaps 10304(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 217
               movddup 5088(pB), rB1
            #endif

         #endif
         #if KB > 217
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10320(pA), rA2
               prefetcht0 5096+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5096(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 218
               movddup 5104(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 218
               movaps 10336(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 218
               movaps 10352(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 218
               movddup 5112(pB), rB1
            #endif

         #endif
         #if KB > 218
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10368(pA), rA2
               prefetcht0 5120(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5120(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 219
               movddup 5128(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 219
               movaps 10384(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 219
               movaps 10400(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 219
               movddup 5136(pB), rB1
            #endif

         #endif
         #if KB > 219
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10416(pA), rA2
               prefetcht0 5144+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5144(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 220
               movddup 5152(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 220
               movaps 10432(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 220
               movaps 10448(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 220
               movddup 5160(pB), rB1
            #endif

         #endif
         #if KB > 220
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10464(pA), rA2
               prefetcht0 5168(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5168(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 221
               movddup 5176(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 221
               movaps 10480(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 221
               movaps 10496(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 221
               movddup 5184(pB), rB1
            #endif

         #endif
         #if KB > 221
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10512(pA), rA2
               prefetcht0 5192+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5192(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 222
               movddup 5200(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 222
               movaps 10528(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 222
               movaps 10544(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 222
               movddup 5208(pB), rB1
            #endif

         #endif
         #if KB > 222
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10560(pA), rA2
               prefetcht0 5216(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5216(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 223
               movddup 5224(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 223
               movaps 10576(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 223
               movaps 10592(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 223
               movddup 5232(pB), rB1
            #endif

         #endif
         #if KB > 223
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10608(pA), rA2
               prefetcht0 5240+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5240(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 224
               movddup 5248(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 224
               movaps 10624(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 224
               movaps 10640(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 224
               movddup 5256(pB), rB1
            #endif

         #endif
         #if KB > 224
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10656(pA), rA2
               prefetcht0 5264(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5264(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 225
               movddup 5272(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 225
               movaps 10672(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 225
               movaps 10688(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 225
               movddup 5280(pB), rB1
            #endif

         #endif
         #if KB > 225
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10704(pA), rA2
               prefetcht0 5288+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5288(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 226
               movddup 5296(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 226
               movaps 10720(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 226
               movaps 10736(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 226
               movddup 5304(pB), rB1
            #endif

         #endif
         #if KB > 226
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10752(pA), rA2
               prefetcht0 5312(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5312(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 227
               movddup 5320(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 227
               movaps 10768(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 227
               movaps 10784(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 227
               movddup 5328(pB), rB1
            #endif

         #endif
         #if KB > 227
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10800(pA), rA2
               prefetcht0 5336+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5336(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 228
               movddup 5344(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 228
               movaps 10816(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 228
               movaps 10832(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 228
               movddup 5352(pB), rB1
            #endif

         #endif
         #if KB > 228
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10848(pA), rA2
               prefetcht0 5360(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5360(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 229
               movddup 5368(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 229
               movaps 10864(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 229
               movaps 10880(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 229
               movddup 5376(pB), rB1
            #endif

         #endif
         #if KB > 229
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10896(pA), rA2
               prefetcht0 5384+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5384(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 230
               movddup 5392(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 230
               movaps 10912(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 230
               movaps 10928(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 230
               movddup 5400(pB), rB1
            #endif

         #endif
         #if KB > 230
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10944(pA), rA2
               prefetcht0 5408(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5408(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 231
               movddup 5416(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 231
               movaps 10960(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 231
               movaps 10976(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 231
               movddup 5424(pB), rB1
            #endif

         #endif
         #if KB > 231
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10992(pA), rA2
               prefetcht0 5432+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5432(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 232
               movddup 5440(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 232
               movaps 11008(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 232
               movaps 11024(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 232
               movddup 5448(pB), rB1
            #endif

         #endif
         #if KB > 232
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11040(pA), rA2
               prefetcht0 5456(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5456(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 233
               movddup 5464(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 233
               movaps 11056(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 233
               movaps 11072(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 233
               movddup 5472(pB), rB1
            #endif

         #endif
         #if KB > 233
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11088(pA), rA2
               prefetcht0 5480+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5480(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 234
               movddup 5488(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 234
               movaps 11104(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 234
               movaps 11120(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 234
               movddup 5496(pB), rB1
            #endif

         #endif
         #if KB > 234
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11136(pA), rA2
               prefetcht0 5504(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5504(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 235
               movddup 5512(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 235
               movaps 11152(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 235
               movaps 11168(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 235
               movddup 5520(pB), rB1
            #endif

         #endif
         #if KB > 235
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11184(pA), rA2
               prefetcht0 5528+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5528(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 236
               movddup 5536(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 236
               movaps 11200(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 236
               movaps 11216(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 236
               movddup 5544(pB), rB1
            #endif

         #endif
         #if KB > 236
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11232(pA), rA2
               prefetcht0 5552(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5552(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 237
               movddup 5560(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 237
               movaps 11248(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 237
               movaps 11264(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 237
               movddup 5568(pB), rB1
            #endif

         #endif
         #if KB > 237
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11280(pA), rA2
               prefetcht0 5576+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5576(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 238
               movddup 5584(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 238
               movaps 11296(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 238
               movaps 11312(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 238
               movddup 5592(pB), rB1
            #endif

         #endif
         #if KB > 238
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11328(pA), rA2
               prefetcht0 5600(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5600(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 239
               movddup 5608(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 239
               movaps 11344(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 239
               movaps 11360(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 239
               movddup 5616(pB), rB1
            #endif

         #endif
         #if KB > 239
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11376(pA), rA2
               prefetcht0 5624+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5624(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 240
               movddup 5632(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 240
               movaps 11392(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 240
               movaps 11408(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 240
               movddup 5640(pB), rB1
            #endif

         #endif
         #if KB > 240
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11424(pA), rA2
               prefetcht0 5648(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5648(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 241
               movddup 5656(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 241
               movaps 11440(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 241
               movaps 11456(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 241
               movddup 5664(pB), rB1
            #endif

         #endif
         #if KB > 241
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11472(pA), rA2
               prefetcht0 5672+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5672(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 242
               movddup 5680(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 242
               movaps 11488(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 242
               movaps 11504(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 242
               movddup 5688(pB), rB1
            #endif

         #endif
         #if KB > 242
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11520(pA), rA2
               prefetcht0 5696(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5696(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 243
               movddup 5704(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 243
               movaps 11536(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 243
               movaps 11552(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 243
               movddup 5712(pB), rB1
            #endif

         #endif
         #if KB > 243
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11568(pA), rA2
               prefetcht0 5720+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5720(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 244
               movddup 5728(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 244
               movaps 11584(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 244
               movaps 11600(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 244
               movddup 5736(pB), rB1
            #endif

         #endif
         #if KB > 244
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11616(pA), rA2
               prefetcht0 5744(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5744(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 245
               movddup 5752(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 245
               movaps 11632(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 245
               movaps 11648(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 245
               movddup 5760(pB), rB1
            #endif

         #endif
         #if KB > 245
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11664(pA), rA2
               prefetcht0 5768+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5768(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 246
               movddup 5776(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 246
               movaps 11680(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 246
               movaps 11696(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 246
               movddup 5784(pB), rB1
            #endif

         #endif
         #if KB > 246
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11712(pA), rA2
               prefetcht0 5792(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5792(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 247
               movddup 5800(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 247
               movaps 11728(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 247
               movaps 11744(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 247
               movddup 5808(pB), rB1
            #endif

         #endif
         #if KB > 247
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11760(pA), rA2
               prefetcht0 5816+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5816(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 248
               movddup 5824(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 248
               movaps 11776(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 248
               movaps 11792(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 248
               movddup 5832(pB), rB1
            #endif

         #endif
         #if KB > 248
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11808(pA), rA2
               prefetcht0 5840(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5840(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 249
               movddup 5848(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 249
               movaps 11824(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 249
               movaps 11840(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 249
               movddup 5856(pB), rB1
            #endif

         #endif
         #if KB > 249
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11856(pA), rA2
               prefetcht0 5864+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5864(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 250
               movddup 5872(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 250
               movaps 11872(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 250
               movaps 11888(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 250
               movddup 5880(pB), rB1
            #endif

         #endif
         #if KB > 250
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11904(pA), rA2
               prefetcht0 5888(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5888(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 251
               movddup 5896(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 251
               movaps 11920(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 251
               movaps 11936(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 251
               movddup 5904(pB), rB1
            #endif

         #endif
         #if KB > 251
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11952(pA), rA2
               prefetcht0 5912+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5912(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 252
               movddup 5920(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 252
               movaps 11968(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 252
               movaps 11984(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 252
               movddup 5928(pB), rB1
            #endif

         #endif
         #if KB > 252
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 12000(pA), rA2
               prefetcht0 5936(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5936(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 253
               movddup 5944(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 253
               movaps 12016(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 253
               movaps 12032(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 253
               movddup 5952(pB), rB1
            #endif

         #endif
         #if KB > 253
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 12048(pA), rA2
               prefetcht0 5960+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5960(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 254
               movddup 5968(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 254
               movaps 12064(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 254
               movaps 12080(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 254
               movddup 5976(pB), rB1
            #endif

         #endif
         #if KB > 254
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 12096(pA), rA2
               prefetcht0 5984(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 5984(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 255
               movddup 5992(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 255
               movaps 12112(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 255
               movaps 12128(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 255
               movddup 6000(pB), rB1
            #endif

         #endif
         #if KB > 255
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 12144(pA), rA2
               prefetcht0 6008+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movddup 6008(pB), rB2
            mulps rA2, rB0
            addps rB0, rC20

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps rA2, rB1
            addps rB1, rC21
            #if KB > 256
               movddup 6016(pB), rB0
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 256
               movaps 12160(pA), rA0
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 256
               movaps 12176(pA), rA1
            #endif

            mulps rB2, rA2
            addps rA2, rC22
            #if KB > 256
               movddup 6024(pB), rB1
            #endif

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
         add $KB*3*8, pB
         add $144, pC
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
