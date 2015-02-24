#include "atlas_asm.h"

#define rm0     %xmm0
#define rB0     %xmm1
#define rB1     %xmm2
#define rB2     %xmm3
#define rA0     %xmm4
#define rA1     %xmm5
#define rA2     %xmm8
#define rC00    %xmm7
#define rC10    %xmm6
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
#if KB <= 128
   #define prefA(m_)
   #define prefB(m_) prefetcht0 m_
#else
   #define prefA(m_)
   #define prefB(m_) prefetcht2 m_
#endif
//#define prefC(m_) prefetcht0 m_
#define prefC(m_) prefetchw m_
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
 * incAm = 6*sizeof*K = 6*8*K = 3*2*8*K = 3*K*16
 */
   lea (KK, KK,2), incAm   /* incAm = 3*K */
   shl $4, incAm           /* incAm = 32*3*K */

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
         movaps -128(pB), rC00
         movaps -128(pA), rA0
         movaps rC00, rC10

         mulps rA0, rC00
         movaps rC10, rC20
         movaps -112(pA), rA1

         mulps rA1, rC10
         movaps -112(pB), rC01
         movaps -96(pA), rA2

         mulps rA2, rC20
         movaps rC01, rC11
         movaps -96(pB), rC02

         mulps rA0, rC01
         movaps rC11, rC21
         movaps -80(pB), rB0

         movaps rC02, rC12
         mulps rA1, rC11
         movaps -64(pB), rB1

         mulps rA2, rC21
         movaps rC12, rC22
         prefC(-128(pC))

         mulps rA0, rC02
         movaps -80(pA), rA0
         prefC(-64(pC))

         mulps rA1, rC12
         movaps -64(pA), rA1
         movaps -48(pB), rB2

         mulps rA2, rC22
         prefC((pC))
         prefC(64(pC))

         #if KB > 1
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps -48(pA), rB0
            addps rB0, rC20
            #if KB > 2
               movaps -32(pB), rB0
            #elif KB == 2
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps -48(pA), rB1
            addps rB1, rC21
            #if KB > 2
               movaps -16(pB), rB1
            #elif KB == 2
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 2
               movaps -32(pA), rA0
            #elif KB == 2
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 2
               movaps -16(pA), rA1
            #elif KB == 2
               prefA(-128(pfA))
            #endif

            mulps -48(pA), rB2
            addps rB2, rC22
            #if KB > 2
               movaps 0(pB), rB2
            #elif KB == 2
               prefA((pfA))
            #endif

         #endif
         #if KB > 2
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 0(pA), rB0
            addps rB0, rC20
            #if KB > 3
               movaps 16(pB), rB0
            #elif KB == 3
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 0(pA), rB1
            addps rB1, rC21
            #if KB > 3
               movaps 32(pB), rB1
            #elif KB == 3
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 3
               movaps 16(pA), rA0
            #elif KB == 3
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 3
               movaps 32(pA), rA1
            #elif KB == 3
               prefA(-128(pfA))
            #endif

            mulps 0(pA), rB2
            addps rB2, rC22
            #if KB > 3
               movaps 48(pB), rB2
            #elif KB == 3
               prefA((pfA))
            #endif

         #endif
         #if KB > 3
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 48(pA), rB0
            addps rB0, rC20
            #if KB > 4
               movaps 64(pB), rB0
            #elif KB == 4
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 48(pA), rB1
            addps rB1, rC21
            #if KB > 4
               movaps 80(pB), rB1
            #elif KB == 4
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 4
               movaps 64(pA), rA0
            #elif KB == 4
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 4
               movaps 80(pA), rA1
            #elif KB == 4
               prefA(-128(pfA))
            #endif

            mulps 48(pA), rB2
            addps rB2, rC22
            #if KB > 4
               movaps 96(pB), rB2
            #elif KB == 4
               prefA((pfA))
            #endif

         #endif
         #if KB > 4
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 96(pA), rB0
            addps rB0, rC20
            #if KB > 5
               movaps 112(pB), rB0
            #elif KB == 5
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 96(pA), rB1
            addps rB1, rC21
            #if KB > 5
               movaps 128(pB), rB1
            #elif KB == 5
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 5
               movaps 112(pA), rA0
            #elif KB == 5
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 5
               movaps 128(pA), rA1
            #elif KB == 5
               prefA(-128(pfA))
            #endif

            mulps 96(pA), rB2
            addps rB2, rC22
            #if KB > 5
               movaps 144(pB), rB2
            #elif KB == 5
               prefA((pfA))
            #endif

         #endif
         #if KB > 5
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 144(pA), rB0
            addps rB0, rC20
            #if KB > 6
               movaps 160(pB), rB0
            #elif KB == 6
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 144(pA), rB1
            addps rB1, rC21
            #if KB > 6
               movaps 176(pB), rB1
            #elif KB == 6
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 6
               movaps 160(pA), rA0
            #elif KB == 6
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 6
               movaps 176(pA), rA1
            #elif KB == 6
               prefA(-128(pfA))
            #endif

            mulps 144(pA), rB2
            addps rB2, rC22
            #if KB > 6
               movaps 192(pB), rB2
            #elif KB == 6
               prefA((pfA))
            #endif

         #endif
         #if KB > 6
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 192(pA), rB0
            addps rB0, rC20
            #if KB > 7
               movaps 208(pB), rB0
            #elif KB == 7
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 192(pA), rB1
            addps rB1, rC21
            #if KB > 7
               movaps 224(pB), rB1
            #elif KB == 7
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 7
               movaps 208(pA), rA0
            #elif KB == 7
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 7
               movaps 224(pA), rA1
            #elif KB == 7
               prefA(-128(pfA))
            #endif

            mulps 192(pA), rB2
            addps rB2, rC22
            #if KB > 7
               movaps 240(pB), rB2
            #elif KB == 7
               prefA((pfA))
            #endif

         #endif
         #if KB > 7
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 240(pA), rB0
            addps rB0, rC20
            #if KB > 8
               movaps 256(pB), rB0
            #elif KB == 8
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 240(pA), rB1
            addps rB1, rC21
            #if KB > 8
               movaps 272(pB), rB1
            #elif KB == 8
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 8
               movaps 256(pA), rA0
            #elif KB == 8
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 8
               movaps 272(pA), rA1
            #elif KB == 8
               prefA(-128(pfA))
            #endif

            mulps 240(pA), rB2
            addps rB2, rC22
            #if KB > 8
               movaps 288(pB), rB2
            #elif KB == 8
               prefA((pfA))
            #endif

         #endif
         #if KB > 8
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 288(pA), rB0
            addps rB0, rC20
            #if KB > 9
               movaps 304(pB), rB0
            #elif KB == 9
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 288(pA), rB1
            addps rB1, rC21
            #if KB > 9
               movaps 320(pB), rB1
            #elif KB == 9
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 9
               movaps 304(pA), rA0
            #elif KB == 9
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 9
               movaps 320(pA), rA1
            #elif KB == 9
               prefA(-128(pfA))
            #endif

            mulps 288(pA), rB2
            addps rB2, rC22
            #if KB > 9
               movaps 336(pB), rB2
            #elif KB == 9
               prefA((pfA))
            #endif

         #endif
         #if KB > 9
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 336(pA), rB0
            addps rB0, rC20
            #if KB > 10
               movaps 352(pB), rB0
            #elif KB == 10
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 336(pA), rB1
            addps rB1, rC21
            #if KB > 10
               movaps 368(pB), rB1
            #elif KB == 10
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 10
               movaps 352(pA), rA0
            #elif KB == 10
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 10
               movaps 368(pA), rA1
            #elif KB == 10
               prefA(-128(pfA))
            #endif

            mulps 336(pA), rB2
            addps rB2, rC22
            #if KB > 10
               movaps 384(pB), rB2
            #elif KB == 10
               prefA((pfA))
            #endif

         #endif
         #if KB > 10
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 384(pA), rB0
            addps rB0, rC20
            #if KB > 11
               movaps 400(pB), rB0
            #elif KB == 11
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 384(pA), rB1
            addps rB1, rC21
            #if KB > 11
               movaps 416(pB), rB1
            #elif KB == 11
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 11
               movaps 400(pA), rA0
            #elif KB == 11
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 11
               movaps 416(pA), rA1
            #elif KB == 11
               prefA(-128(pfA))
            #endif

            mulps 384(pA), rB2
            addps rB2, rC22
            #if KB > 11
               movaps 432(pB), rB2
            #elif KB == 11
               prefA((pfA))
            #endif

         #endif
         #if KB > 11
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 432(pA), rB0
            addps rB0, rC20
            #if KB > 12
               movaps 448(pB), rB0
            #elif KB == 12
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 432(pA), rB1
            addps rB1, rC21
            #if KB > 12
               movaps 464(pB), rB1
            #elif KB == 12
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 12
               movaps 448(pA), rA0
            #elif KB == 12
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 12
               movaps 464(pA), rA1
            #elif KB == 12
               prefA(-128(pfA))
            #endif

            mulps 432(pA), rB2
            addps rB2, rC22
            #if KB > 12
               movaps 480(pB), rB2
            #elif KB == 12
               prefA((pfA))
            #endif

         #endif
         #if KB > 12
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 480(pA), rB0
            addps rB0, rC20
            #if KB > 13
               movaps 496(pB), rB0
            #elif KB == 13
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 480(pA), rB1
            addps rB1, rC21
            #if KB > 13
               movaps 512(pB), rB1
            #elif KB == 13
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 13
               movaps 496(pA), rA0
            #elif KB == 13
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 13
               movaps 512(pA), rA1
            #elif KB == 13
               prefA(-128(pfA))
            #endif

            mulps 480(pA), rB2
            addps rB2, rC22
            #if KB > 13
               movaps 528(pB), rB2
            #elif KB == 13
               prefA((pfA))
            #endif

         #endif
         #if KB > 13
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 528(pA), rB0
            addps rB0, rC20
            #if KB > 14
               movaps 544(pB), rB0
            #elif KB == 14
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 528(pA), rB1
            addps rB1, rC21
            #if KB > 14
               movaps 560(pB), rB1
            #elif KB == 14
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 14
               movaps 544(pA), rA0
            #elif KB == 14
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 14
               movaps 560(pA), rA1
            #elif KB == 14
               prefA(-128(pfA))
            #endif

            mulps 528(pA), rB2
            addps rB2, rC22
            #if KB > 14
               movaps 576(pB), rB2
            #elif KB == 14
               prefA((pfA))
            #endif

         #endif
         #if KB > 14
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 576(pA), rB0
            addps rB0, rC20
            #if KB > 15
               movaps 592(pB), rB0
            #elif KB == 15
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 576(pA), rB1
            addps rB1, rC21
            #if KB > 15
               movaps 608(pB), rB1
            #elif KB == 15
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 15
               movaps 592(pA), rA0
            #elif KB == 15
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 15
               movaps 608(pA), rA1
            #elif KB == 15
               prefA(-128(pfA))
            #endif

            mulps 576(pA), rB2
            addps rB2, rC22
            #if KB > 15
               movaps 624(pB), rB2
            #elif KB == 15
               prefA((pfA))
            #endif

         #endif
         #if KB > 15
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 624(pA), rB0
            addps rB0, rC20
            #if KB > 16
               movaps 640(pB), rB0
            #elif KB == 16
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 624(pA), rB1
            addps rB1, rC21
            #if KB > 16
               movaps 656(pB), rB1
            #elif KB == 16
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 16
               movaps 640(pA), rA0
            #elif KB == 16
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 16
               movaps 656(pA), rA1
            #elif KB == 16
               prefA(-128(pfA))
            #endif

            mulps 624(pA), rB2
            addps rB2, rC22
            #if KB > 16
               movaps 672(pB), rB2
            #elif KB == 16
               prefA((pfA))
            #endif

         #endif
         #if KB > 16
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 672(pA), rB0
            addps rB0, rC20
            #if KB > 17
               movaps 688(pB), rB0
            #elif KB == 17
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 672(pA), rB1
            addps rB1, rC21
            #if KB > 17
               movaps 704(pB), rB1
            #elif KB == 17
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 17
               movaps 688(pA), rA0
            #elif KB == 17
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 17
               movaps 704(pA), rA1
            #elif KB == 17
               prefA(-128(pfA))
            #endif

            mulps 672(pA), rB2
            addps rB2, rC22
            #if KB > 17
               movaps 720(pB), rB2
            #elif KB == 17
               prefA((pfA))
            #endif

         #endif
         #if KB > 17
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 720(pA), rB0
            addps rB0, rC20
            #if KB > 18
               movaps 736(pB), rB0
            #elif KB == 18
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 720(pA), rB1
            addps rB1, rC21
            #if KB > 18
               movaps 752(pB), rB1
            #elif KB == 18
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 18
               movaps 736(pA), rA0
            #elif KB == 18
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 18
               movaps 752(pA), rA1
            #elif KB == 18
               prefA(-128(pfA))
            #endif

            mulps 720(pA), rB2
            addps rB2, rC22
            #if KB > 18
               movaps 768(pB), rB2
            #elif KB == 18
               prefA((pfA))
            #endif

         #endif
         #if KB > 18
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 768(pA), rB0
            addps rB0, rC20
            #if KB > 19
               movaps 784(pB), rB0
            #elif KB == 19
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 768(pA), rB1
            addps rB1, rC21
            #if KB > 19
               movaps 800(pB), rB1
            #elif KB == 19
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 19
               movaps 784(pA), rA0
            #elif KB == 19
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 19
               movaps 800(pA), rA1
            #elif KB == 19
               prefA(-128(pfA))
            #endif

            mulps 768(pA), rB2
            addps rB2, rC22
            #if KB > 19
               movaps 816(pB), rB2
            #elif KB == 19
               prefA((pfA))
            #endif

         #endif
         #if KB > 19
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 816(pA), rB0
            addps rB0, rC20
            #if KB > 20
               movaps 832(pB), rB0
            #elif KB == 20
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 816(pA), rB1
            addps rB1, rC21
            #if KB > 20
               movaps 848(pB), rB1
            #elif KB == 20
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 20
               movaps 832(pA), rA0
            #elif KB == 20
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 20
               movaps 848(pA), rA1
            #elif KB == 20
               prefA(-128(pfA))
            #endif

            mulps 816(pA), rB2
            addps rB2, rC22
            #if KB > 20
               movaps 864(pB), rB2
            #elif KB == 20
               prefA((pfA))
            #endif

         #endif
         #if KB > 20
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 864(pA), rB0
            addps rB0, rC20
            #if KB > 21
               movaps 880(pB), rB0
            #elif KB == 21
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 864(pA), rB1
            addps rB1, rC21
            #if KB > 21
               movaps 896(pB), rB1
            #elif KB == 21
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 21
               movaps 880(pA), rA0
            #elif KB == 21
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 21
               movaps 896(pA), rA1
            #elif KB == 21
               prefA(-128(pfA))
            #endif

            mulps 864(pA), rB2
            addps rB2, rC22
            #if KB > 21
               movaps 912(pB), rB2
            #elif KB == 21
               prefA((pfA))
            #endif

         #endif
         #if KB > 21
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 912(pA), rB0
            addps rB0, rC20
            #if KB > 22
               movaps 928(pB), rB0
            #elif KB == 22
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 912(pA), rB1
            addps rB1, rC21
            #if KB > 22
               movaps 944(pB), rB1
            #elif KB == 22
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 22
               movaps 928(pA), rA0
            #elif KB == 22
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 22
               movaps 944(pA), rA1
            #elif KB == 22
               prefA(-128(pfA))
            #endif

            mulps 912(pA), rB2
            addps rB2, rC22
            #if KB > 22
               movaps 960(pB), rB2
            #elif KB == 22
               prefA((pfA))
            #endif

         #endif
         #if KB > 22
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 960(pA), rB0
            addps rB0, rC20
            #if KB > 23
               movaps 976(pB), rB0
            #elif KB == 23
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 960(pA), rB1
            addps rB1, rC21
            #if KB > 23
               movaps 992(pB), rB1
            #elif KB == 23
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 23
               movaps 976(pA), rA0
            #elif KB == 23
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 23
               movaps 992(pA), rA1
            #elif KB == 23
               prefA(-128(pfA))
            #endif

            mulps 960(pA), rB2
            addps rB2, rC22
            #if KB > 23
               movaps 1008(pB), rB2
            #elif KB == 23
               prefA((pfA))
            #endif

         #endif
         #if KB > 23
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1008(pA), rB0
            addps rB0, rC20
            #if KB > 24
               movaps 1024(pB), rB0
            #elif KB == 24
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1008(pA), rB1
            addps rB1, rC21
            #if KB > 24
               movaps 1040(pB), rB1
            #elif KB == 24
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 24
               movaps 1024(pA), rA0
            #elif KB == 24
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 24
               movaps 1040(pA), rA1
            #elif KB == 24
               prefA(-128(pfA))
            #endif

            mulps 1008(pA), rB2
            addps rB2, rC22
            #if KB > 24
               movaps 1056(pB), rB2
            #elif KB == 24
               prefA((pfA))
            #endif

         #endif
         #if KB > 24
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1056(pA), rB0
            addps rB0, rC20
            #if KB > 25
               movaps 1072(pB), rB0
            #elif KB == 25
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1056(pA), rB1
            addps rB1, rC21
            #if KB > 25
               movaps 1088(pB), rB1
            #elif KB == 25
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 25
               movaps 1072(pA), rA0
            #elif KB == 25
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 25
               movaps 1088(pA), rA1
            #elif KB == 25
               prefA(-128(pfA))
            #endif

            mulps 1056(pA), rB2
            addps rB2, rC22
            #if KB > 25
               movaps 1104(pB), rB2
            #elif KB == 25
               prefA((pfA))
            #endif

         #endif
         #if KB > 25
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1104(pA), rB0
            addps rB0, rC20
            #if KB > 26
               movaps 1120(pB), rB0
            #elif KB == 26
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1104(pA), rB1
            addps rB1, rC21
            #if KB > 26
               movaps 1136(pB), rB1
            #elif KB == 26
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 26
               movaps 1120(pA), rA0
            #elif KB == 26
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 26
               movaps 1136(pA), rA1
            #elif KB == 26
               prefA(-128(pfA))
            #endif

            mulps 1104(pA), rB2
            addps rB2, rC22
            #if KB > 26
               movaps 1152(pB), rB2
            #elif KB == 26
               prefA((pfA))
            #endif

         #endif
         #if KB > 26
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1152(pA), rB0
            addps rB0, rC20
            #if KB > 27
               movaps 1168(pB), rB0
            #elif KB == 27
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1152(pA), rB1
            addps rB1, rC21
            #if KB > 27
               movaps 1184(pB), rB1
            #elif KB == 27
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 27
               movaps 1168(pA), rA0
            #elif KB == 27
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 27
               movaps 1184(pA), rA1
            #elif KB == 27
               prefA(-128(pfA))
            #endif

            mulps 1152(pA), rB2
            addps rB2, rC22
            #if KB > 27
               movaps 1200(pB), rB2
            #elif KB == 27
               prefA((pfA))
            #endif

         #endif
         #if KB > 27
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1200(pA), rB0
            addps rB0, rC20
            #if KB > 28
               movaps 1216(pB), rB0
            #elif KB == 28
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1200(pA), rB1
            addps rB1, rC21
            #if KB > 28
               movaps 1232(pB), rB1
            #elif KB == 28
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 28
               movaps 1216(pA), rA0
            #elif KB == 28
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 28
               movaps 1232(pA), rA1
            #elif KB == 28
               prefA(-128(pfA))
            #endif

            mulps 1200(pA), rB2
            addps rB2, rC22
            #if KB > 28
               movaps 1248(pB), rB2
            #elif KB == 28
               prefA((pfA))
            #endif

         #endif
         #if KB > 28
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1248(pA), rB0
            addps rB0, rC20
            #if KB > 29
               movaps 1264(pB), rB0
            #elif KB == 29
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1248(pA), rB1
            addps rB1, rC21
            #if KB > 29
               movaps 1280(pB), rB1
            #elif KB == 29
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 29
               movaps 1264(pA), rA0
            #elif KB == 29
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 29
               movaps 1280(pA), rA1
            #elif KB == 29
               prefA(-128(pfA))
            #endif

            mulps 1248(pA), rB2
            addps rB2, rC22
            #if KB > 29
               movaps 1296(pB), rB2
            #elif KB == 29
               prefA((pfA))
            #endif

         #endif
         #if KB > 29
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1296(pA), rB0
            addps rB0, rC20
            #if KB > 30
               movaps 1312(pB), rB0
            #elif KB == 30
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1296(pA), rB1
            addps rB1, rC21
            #if KB > 30
               movaps 1328(pB), rB1
            #elif KB == 30
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 30
               movaps 1312(pA), rA0
            #elif KB == 30
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 30
               movaps 1328(pA), rA1
            #elif KB == 30
               prefA(-128(pfA))
            #endif

            mulps 1296(pA), rB2
            addps rB2, rC22
            #if KB > 30
               movaps 1344(pB), rB2
            #elif KB == 30
               prefA((pfA))
            #endif

         #endif
         #if KB > 30
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1344(pA), rB0
            addps rB0, rC20
            #if KB > 31
               movaps 1360(pB), rB0
            #elif KB == 31
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1344(pA), rB1
            addps rB1, rC21
            #if KB > 31
               movaps 1376(pB), rB1
            #elif KB == 31
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 31
               movaps 1360(pA), rA0
            #elif KB == 31
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 31
               movaps 1376(pA), rA1
            #elif KB == 31
               prefA(-128(pfA))
            #endif

            mulps 1344(pA), rB2
            addps rB2, rC22
            #if KB > 31
               movaps 1392(pB), rB2
            #elif KB == 31
               prefA((pfA))
            #endif

         #endif
         #if KB > 31
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1392(pA), rB0
            addps rB0, rC20
            #if KB > 32
               movaps 1408(pB), rB0
            #elif KB == 32
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1392(pA), rB1
            addps rB1, rC21
            #if KB > 32
               movaps 1424(pB), rB1
            #elif KB == 32
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 32
               movaps 1408(pA), rA0
            #elif KB == 32
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 32
               movaps 1424(pA), rA1
            #elif KB == 32
               prefA(-128(pfA))
            #endif

            mulps 1392(pA), rB2
            addps rB2, rC22
            #if KB > 32
               movaps 1440(pB), rB2
            #elif KB == 32
               prefA((pfA))
            #endif

         #endif
         #if KB > 32
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1440(pA), rB0
            addps rB0, rC20
            #if KB > 33
               movaps 1456(pB), rB0
            #elif KB == 33
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1440(pA), rB1
            addps rB1, rC21
            #if KB > 33
               movaps 1472(pB), rB1
            #elif KB == 33
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 33
               movaps 1456(pA), rA0
            #elif KB == 33
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 33
               movaps 1472(pA), rA1
            #elif KB == 33
               prefA(-128(pfA))
            #endif

            mulps 1440(pA), rB2
            addps rB2, rC22
            #if KB > 33
               movaps 1488(pB), rB2
            #elif KB == 33
               prefA((pfA))
            #endif

         #endif
         #if KB > 33
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1488(pA), rB0
            addps rB0, rC20
            #if KB > 34
               movaps 1504(pB), rB0
            #elif KB == 34
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1488(pA), rB1
            addps rB1, rC21
            #if KB > 34
               movaps 1520(pB), rB1
            #elif KB == 34
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 34
               movaps 1504(pA), rA0
            #elif KB == 34
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 34
               movaps 1520(pA), rA1
            #elif KB == 34
               prefA(-128(pfA))
            #endif

            mulps 1488(pA), rB2
            addps rB2, rC22
            #if KB > 34
               movaps 1536(pB), rB2
            #elif KB == 34
               prefA((pfA))
            #endif

         #endif
         #if KB > 34
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1536(pA), rB0
            addps rB0, rC20
            #if KB > 35
               movaps 1552(pB), rB0
            #elif KB == 35
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1536(pA), rB1
            addps rB1, rC21
            #if KB > 35
               movaps 1568(pB), rB1
            #elif KB == 35
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 35
               movaps 1552(pA), rA0
            #elif KB == 35
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 35
               movaps 1568(pA), rA1
            #elif KB == 35
               prefA(-128(pfA))
            #endif

            mulps 1536(pA), rB2
            addps rB2, rC22
            #if KB > 35
               movaps 1584(pB), rB2
            #elif KB == 35
               prefA((pfA))
            #endif

         #endif
         #if KB > 35
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1584(pA), rB0
            addps rB0, rC20
            #if KB > 36
               movaps 1600(pB), rB0
            #elif KB == 36
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1584(pA), rB1
            addps rB1, rC21
            #if KB > 36
               movaps 1616(pB), rB1
            #elif KB == 36
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 36
               movaps 1600(pA), rA0
            #elif KB == 36
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 36
               movaps 1616(pA), rA1
            #elif KB == 36
               prefA(-128(pfA))
            #endif

            mulps 1584(pA), rB2
            addps rB2, rC22
            #if KB > 36
               movaps 1632(pB), rB2
            #elif KB == 36
               prefA((pfA))
            #endif

         #endif
         #if KB > 36
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1632(pA), rB0
            addps rB0, rC20
            #if KB > 37
               movaps 1648(pB), rB0
            #elif KB == 37
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1632(pA), rB1
            addps rB1, rC21
            #if KB > 37
               movaps 1664(pB), rB1
            #elif KB == 37
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 37
               movaps 1648(pA), rA0
            #elif KB == 37
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 37
               movaps 1664(pA), rA1
            #elif KB == 37
               prefA(-128(pfA))
            #endif

            mulps 1632(pA), rB2
            addps rB2, rC22
            #if KB > 37
               movaps 1680(pB), rB2
            #elif KB == 37
               prefA((pfA))
            #endif

         #endif
         #if KB > 37
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1680(pA), rB0
            addps rB0, rC20
            #if KB > 38
               movaps 1696(pB), rB0
            #elif KB == 38
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1680(pA), rB1
            addps rB1, rC21
            #if KB > 38
               movaps 1712(pB), rB1
            #elif KB == 38
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 38
               movaps 1696(pA), rA0
            #elif KB == 38
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 38
               movaps 1712(pA), rA1
            #elif KB == 38
               prefA(-128(pfA))
            #endif

            mulps 1680(pA), rB2
            addps rB2, rC22
            #if KB > 38
               movaps 1728(pB), rB2
            #elif KB == 38
               prefA((pfA))
            #endif

         #endif
         #if KB > 38
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1728(pA), rB0
            addps rB0, rC20
            #if KB > 39
               movaps 1744(pB), rB0
            #elif KB == 39
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1728(pA), rB1
            addps rB1, rC21
            #if KB > 39
               movaps 1760(pB), rB1
            #elif KB == 39
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 39
               movaps 1744(pA), rA0
            #elif KB == 39
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 39
               movaps 1760(pA), rA1
            #elif KB == 39
               prefA(-128(pfA))
            #endif

            mulps 1728(pA), rB2
            addps rB2, rC22
            #if KB > 39
               movaps 1776(pB), rB2
            #elif KB == 39
               prefA((pfA))
            #endif

         #endif
         #if KB > 39
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1776(pA), rB0
            addps rB0, rC20
            #if KB > 40
               movaps 1792(pB), rB0
            #elif KB == 40
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1776(pA), rB1
            addps rB1, rC21
            #if KB > 40
               movaps 1808(pB), rB1
            #elif KB == 40
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 40
               movaps 1792(pA), rA0
            #elif KB == 40
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 40
               movaps 1808(pA), rA1
            #elif KB == 40
               prefA(-128(pfA))
            #endif

            mulps 1776(pA), rB2
            addps rB2, rC22
            #if KB > 40
               movaps 1824(pB), rB2
            #elif KB == 40
               prefA((pfA))
            #endif

         #endif
         #if KB > 40
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1824(pA), rB0
            addps rB0, rC20
            #if KB > 41
               movaps 1840(pB), rB0
            #elif KB == 41
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1824(pA), rB1
            addps rB1, rC21
            #if KB > 41
               movaps 1856(pB), rB1
            #elif KB == 41
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 41
               movaps 1840(pA), rA0
            #elif KB == 41
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 41
               movaps 1856(pA), rA1
            #elif KB == 41
               prefA(-128(pfA))
            #endif

            mulps 1824(pA), rB2
            addps rB2, rC22
            #if KB > 41
               movaps 1872(pB), rB2
            #elif KB == 41
               prefA((pfA))
            #endif

         #endif
         #if KB > 41
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1872(pA), rB0
            addps rB0, rC20
            #if KB > 42
               movaps 1888(pB), rB0
            #elif KB == 42
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1872(pA), rB1
            addps rB1, rC21
            #if KB > 42
               movaps 1904(pB), rB1
            #elif KB == 42
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 42
               movaps 1888(pA), rA0
            #elif KB == 42
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 42
               movaps 1904(pA), rA1
            #elif KB == 42
               prefA(-128(pfA))
            #endif

            mulps 1872(pA), rB2
            addps rB2, rC22
            #if KB > 42
               movaps 1920(pB), rB2
            #elif KB == 42
               prefA((pfA))
            #endif

         #endif
         #if KB > 42
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1920(pA), rB0
            addps rB0, rC20
            #if KB > 43
               movaps 1936(pB), rB0
            #elif KB == 43
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1920(pA), rB1
            addps rB1, rC21
            #if KB > 43
               movaps 1952(pB), rB1
            #elif KB == 43
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 43
               movaps 1936(pA), rA0
            #elif KB == 43
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 43
               movaps 1952(pA), rA1
            #elif KB == 43
               prefA(-128(pfA))
            #endif

            mulps 1920(pA), rB2
            addps rB2, rC22
            #if KB > 43
               movaps 1968(pB), rB2
            #elif KB == 43
               prefA((pfA))
            #endif

         #endif
         #if KB > 43
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 1968(pA), rB0
            addps rB0, rC20
            #if KB > 44
               movaps 1984(pB), rB0
            #elif KB == 44
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 1968(pA), rB1
            addps rB1, rC21
            #if KB > 44
               movaps 2000(pB), rB1
            #elif KB == 44
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 44
               movaps 1984(pA), rA0
            #elif KB == 44
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 44
               movaps 2000(pA), rA1
            #elif KB == 44
               prefA(-128(pfA))
            #endif

            mulps 1968(pA), rB2
            addps rB2, rC22
            #if KB > 44
               movaps 2016(pB), rB2
            #elif KB == 44
               prefA((pfA))
            #endif

         #endif
         #if KB > 44
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2016(pA), rB0
            addps rB0, rC20
            #if KB > 45
               movaps 2032(pB), rB0
            #elif KB == 45
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2016(pA), rB1
            addps rB1, rC21
            #if KB > 45
               movaps 2048(pB), rB1
            #elif KB == 45
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 45
               movaps 2032(pA), rA0
            #elif KB == 45
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 45
               movaps 2048(pA), rA1
            #elif KB == 45
               prefA(-128(pfA))
            #endif

            mulps 2016(pA), rB2
            addps rB2, rC22
            #if KB > 45
               movaps 2064(pB), rB2
            #elif KB == 45
               prefA((pfA))
            #endif

         #endif
         #if KB > 45
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2064(pA), rB0
            addps rB0, rC20
            #if KB > 46
               movaps 2080(pB), rB0
            #elif KB == 46
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2064(pA), rB1
            addps rB1, rC21
            #if KB > 46
               movaps 2096(pB), rB1
            #elif KB == 46
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 46
               movaps 2080(pA), rA0
            #elif KB == 46
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 46
               movaps 2096(pA), rA1
            #elif KB == 46
               prefA(-128(pfA))
            #endif

            mulps 2064(pA), rB2
            addps rB2, rC22
            #if KB > 46
               movaps 2112(pB), rB2
            #elif KB == 46
               prefA((pfA))
            #endif

         #endif
         #if KB > 46
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2112(pA), rB0
            addps rB0, rC20
            #if KB > 47
               movaps 2128(pB), rB0
            #elif KB == 47
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2112(pA), rB1
            addps rB1, rC21
            #if KB > 47
               movaps 2144(pB), rB1
            #elif KB == 47
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 47
               movaps 2128(pA), rA0
            #elif KB == 47
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 47
               movaps 2144(pA), rA1
            #elif KB == 47
               prefA(-128(pfA))
            #endif

            mulps 2112(pA), rB2
            addps rB2, rC22
            #if KB > 47
               movaps 2160(pB), rB2
            #elif KB == 47
               prefA((pfA))
            #endif

         #endif
         #if KB > 47
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2160(pA), rB0
            addps rB0, rC20
            #if KB > 48
               movaps 2176(pB), rB0
            #elif KB == 48
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2160(pA), rB1
            addps rB1, rC21
            #if KB > 48
               movaps 2192(pB), rB1
            #elif KB == 48
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 48
               movaps 2176(pA), rA0
            #elif KB == 48
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 48
               movaps 2192(pA), rA1
            #elif KB == 48
               prefA(-128(pfA))
            #endif

            mulps 2160(pA), rB2
            addps rB2, rC22
            #if KB > 48
               movaps 2208(pB), rB2
            #elif KB == 48
               prefA((pfA))
            #endif

         #endif
         #if KB > 48
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2208(pA), rB0
            addps rB0, rC20
            #if KB > 49
               movaps 2224(pB), rB0
            #elif KB == 49
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2208(pA), rB1
            addps rB1, rC21
            #if KB > 49
               movaps 2240(pB), rB1
            #elif KB == 49
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 49
               movaps 2224(pA), rA0
            #elif KB == 49
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 49
               movaps 2240(pA), rA1
            #elif KB == 49
               prefA(-128(pfA))
            #endif

            mulps 2208(pA), rB2
            addps rB2, rC22
            #if KB > 49
               movaps 2256(pB), rB2
            #elif KB == 49
               prefA((pfA))
            #endif

         #endif
         #if KB > 49
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2256(pA), rB0
            addps rB0, rC20
            #if KB > 50
               movaps 2272(pB), rB0
            #elif KB == 50
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2256(pA), rB1
            addps rB1, rC21
            #if KB > 50
               movaps 2288(pB), rB1
            #elif KB == 50
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 50
               movaps 2272(pA), rA0
            #elif KB == 50
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 50
               movaps 2288(pA), rA1
            #elif KB == 50
               prefA(-128(pfA))
            #endif

            mulps 2256(pA), rB2
            addps rB2, rC22
            #if KB > 50
               movaps 2304(pB), rB2
            #elif KB == 50
               prefA((pfA))
            #endif

         #endif
         #if KB > 50
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2304(pA), rB0
            addps rB0, rC20
            #if KB > 51
               movaps 2320(pB), rB0
            #elif KB == 51
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2304(pA), rB1
            addps rB1, rC21
            #if KB > 51
               movaps 2336(pB), rB1
            #elif KB == 51
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 51
               movaps 2320(pA), rA0
            #elif KB == 51
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 51
               movaps 2336(pA), rA1
            #elif KB == 51
               prefA(-128(pfA))
            #endif

            mulps 2304(pA), rB2
            addps rB2, rC22
            #if KB > 51
               movaps 2352(pB), rB2
            #elif KB == 51
               prefA((pfA))
            #endif

         #endif
         #if KB > 51
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2352(pA), rB0
            addps rB0, rC20
            #if KB > 52
               movaps 2368(pB), rB0
            #elif KB == 52
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2352(pA), rB1
            addps rB1, rC21
            #if KB > 52
               movaps 2384(pB), rB1
            #elif KB == 52
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 52
               movaps 2368(pA), rA0
            #elif KB == 52
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 52
               movaps 2384(pA), rA1
            #elif KB == 52
               prefA(-128(pfA))
            #endif

            mulps 2352(pA), rB2
            addps rB2, rC22
            #if KB > 52
               movaps 2400(pB), rB2
            #elif KB == 52
               prefA((pfA))
            #endif

         #endif
         #if KB > 52
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2400(pA), rB0
            addps rB0, rC20
            #if KB > 53
               movaps 2416(pB), rB0
            #elif KB == 53
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2400(pA), rB1
            addps rB1, rC21
            #if KB > 53
               movaps 2432(pB), rB1
            #elif KB == 53
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 53
               movaps 2416(pA), rA0
            #elif KB == 53
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 53
               movaps 2432(pA), rA1
            #elif KB == 53
               prefA(-128(pfA))
            #endif

            mulps 2400(pA), rB2
            addps rB2, rC22
            #if KB > 53
               movaps 2448(pB), rB2
            #elif KB == 53
               prefA((pfA))
            #endif

         #endif
         #if KB > 53
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2448(pA), rB0
            addps rB0, rC20
            #if KB > 54
               movaps 2464(pB), rB0
            #elif KB == 54
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2448(pA), rB1
            addps rB1, rC21
            #if KB > 54
               movaps 2480(pB), rB1
            #elif KB == 54
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 54
               movaps 2464(pA), rA0
            #elif KB == 54
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 54
               movaps 2480(pA), rA1
            #elif KB == 54
               prefA(-128(pfA))
            #endif

            mulps 2448(pA), rB2
            addps rB2, rC22
            #if KB > 54
               movaps 2496(pB), rB2
            #elif KB == 54
               prefA((pfA))
            #endif

         #endif
         #if KB > 54
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2496(pA), rB0
            addps rB0, rC20
            #if KB > 55
               movaps 2512(pB), rB0
            #elif KB == 55
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2496(pA), rB1
            addps rB1, rC21
            #if KB > 55
               movaps 2528(pB), rB1
            #elif KB == 55
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 55
               movaps 2512(pA), rA0
            #elif KB == 55
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 55
               movaps 2528(pA), rA1
            #elif KB == 55
               prefA(-128(pfA))
            #endif

            mulps 2496(pA), rB2
            addps rB2, rC22
            #if KB > 55
               movaps 2544(pB), rB2
            #elif KB == 55
               prefA((pfA))
            #endif

         #endif
         #if KB > 55
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2544(pA), rB0
            addps rB0, rC20
            #if KB > 56
               movaps 2560(pB), rB0
            #elif KB == 56
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2544(pA), rB1
            addps rB1, rC21
            #if KB > 56
               movaps 2576(pB), rB1
            #elif KB == 56
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 56
               movaps 2560(pA), rA0
            #elif KB == 56
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 56
               movaps 2576(pA), rA1
            #elif KB == 56
               prefA(-128(pfA))
            #endif

            mulps 2544(pA), rB2
            addps rB2, rC22
            #if KB > 56
               movaps 2592(pB), rB2
            #elif KB == 56
               prefA((pfA))
            #endif

         #endif
         #if KB > 56
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2592(pA), rB0
            addps rB0, rC20
            #if KB > 57
               movaps 2608(pB), rB0
            #elif KB == 57
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2592(pA), rB1
            addps rB1, rC21
            #if KB > 57
               movaps 2624(pB), rB1
            #elif KB == 57
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 57
               movaps 2608(pA), rA0
            #elif KB == 57
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 57
               movaps 2624(pA), rA1
            #elif KB == 57
               prefA(-128(pfA))
            #endif

            mulps 2592(pA), rB2
            addps rB2, rC22
            #if KB > 57
               movaps 2640(pB), rB2
            #elif KB == 57
               prefA((pfA))
            #endif

         #endif
         #if KB > 57
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2640(pA), rB0
            addps rB0, rC20
            #if KB > 58
               movaps 2656(pB), rB0
            #elif KB == 58
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2640(pA), rB1
            addps rB1, rC21
            #if KB > 58
               movaps 2672(pB), rB1
            #elif KB == 58
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 58
               movaps 2656(pA), rA0
            #elif KB == 58
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 58
               movaps 2672(pA), rA1
            #elif KB == 58
               prefA(-128(pfA))
            #endif

            mulps 2640(pA), rB2
            addps rB2, rC22
            #if KB > 58
               movaps 2688(pB), rB2
            #elif KB == 58
               prefA((pfA))
            #endif

         #endif
         #if KB > 58
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2688(pA), rB0
            addps rB0, rC20
            #if KB > 59
               movaps 2704(pB), rB0
            #elif KB == 59
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2688(pA), rB1
            addps rB1, rC21
            #if KB > 59
               movaps 2720(pB), rB1
            #elif KB == 59
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 59
               movaps 2704(pA), rA0
            #elif KB == 59
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 59
               movaps 2720(pA), rA1
            #elif KB == 59
               prefA(-128(pfA))
            #endif

            mulps 2688(pA), rB2
            addps rB2, rC22
            #if KB > 59
               movaps 2736(pB), rB2
            #elif KB == 59
               prefA((pfA))
            #endif

         #endif
         #if KB > 59
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2736(pA), rB0
            addps rB0, rC20
            #if KB > 60
               movaps 2752(pB), rB0
            #elif KB == 60
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2736(pA), rB1
            addps rB1, rC21
            #if KB > 60
               movaps 2768(pB), rB1
            #elif KB == 60
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 60
               movaps 2752(pA), rA0
            #elif KB == 60
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 60
               movaps 2768(pA), rA1
            #elif KB == 60
               prefA(-128(pfA))
            #endif

            mulps 2736(pA), rB2
            addps rB2, rC22
            #if KB > 60
               movaps 2784(pB), rB2
            #elif KB == 60
               prefA((pfA))
            #endif

         #endif
         #if KB > 60
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2784(pA), rB0
            addps rB0, rC20
            #if KB > 61
               movaps 2800(pB), rB0
            #elif KB == 61
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2784(pA), rB1
            addps rB1, rC21
            #if KB > 61
               movaps 2816(pB), rB1
            #elif KB == 61
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 61
               movaps 2800(pA), rA0
            #elif KB == 61
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 61
               movaps 2816(pA), rA1
            #elif KB == 61
               prefA(-128(pfA))
            #endif

            mulps 2784(pA), rB2
            addps rB2, rC22
            #if KB > 61
               movaps 2832(pB), rB2
            #elif KB == 61
               prefA((pfA))
            #endif

         #endif
         #if KB > 61
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2832(pA), rB0
            addps rB0, rC20
            #if KB > 62
               movaps 2848(pB), rB0
            #elif KB == 62
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2832(pA), rB1
            addps rB1, rC21
            #if KB > 62
               movaps 2864(pB), rB1
            #elif KB == 62
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 62
               movaps 2848(pA), rA0
            #elif KB == 62
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 62
               movaps 2864(pA), rA1
            #elif KB == 62
               prefA(-128(pfA))
            #endif

            mulps 2832(pA), rB2
            addps rB2, rC22
            #if KB > 62
               movaps 2880(pB), rB2
            #elif KB == 62
               prefA((pfA))
            #endif

         #endif
         #if KB > 62
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2880(pA), rB0
            addps rB0, rC20
            #if KB > 63
               movaps 2896(pB), rB0
            #elif KB == 63
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2880(pA), rB1
            addps rB1, rC21
            #if KB > 63
               movaps 2912(pB), rB1
            #elif KB == 63
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 63
               movaps 2896(pA), rA0
            #elif KB == 63
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 63
               movaps 2912(pA), rA1
            #elif KB == 63
               prefA(-128(pfA))
            #endif

            mulps 2880(pA), rB2
            addps rB2, rC22
            #if KB > 63
               movaps 2928(pB), rB2
            #elif KB == 63
               prefA((pfA))
            #endif

         #endif
         #if KB > 63
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2928(pA), rB0
            addps rB0, rC20
            #if KB > 64
               movaps 2944(pB), rB0
            #elif KB == 64
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2928(pA), rB1
            addps rB1, rC21
            #if KB > 64
               movaps 2960(pB), rB1
            #elif KB == 64
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 64
               movaps 2944(pA), rA0
            #elif KB == 64
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 64
               movaps 2960(pA), rA1
            #elif KB == 64
               prefA(-128(pfA))
            #endif

            mulps 2928(pA), rB2
            addps rB2, rC22
            #if KB > 64
               movaps 2976(pB), rB2
            #elif KB == 64
               prefA((pfA))
            #endif

         #endif
         #if KB > 64
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 2976(pA), rB0
            addps rB0, rC20
            #if KB > 65
               movaps 2992(pB), rB0
            #elif KB == 65
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 2976(pA), rB1
            addps rB1, rC21
            #if KB > 65
               movaps 3008(pB), rB1
            #elif KB == 65
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 65
               movaps 2992(pA), rA0
            #elif KB == 65
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 65
               movaps 3008(pA), rA1
            #elif KB == 65
               prefA(-128(pfA))
            #endif

            mulps 2976(pA), rB2
            addps rB2, rC22
            #if KB > 65
               movaps 3024(pB), rB2
            #elif KB == 65
               prefA((pfA))
            #endif

         #endif
         #if KB > 65
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3024(pA), rB0
            addps rB0, rC20
            #if KB > 66
               movaps 3040(pB), rB0
            #elif KB == 66
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3024(pA), rB1
            addps rB1, rC21
            #if KB > 66
               movaps 3056(pB), rB1
            #elif KB == 66
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 66
               movaps 3040(pA), rA0
            #elif KB == 66
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 66
               movaps 3056(pA), rA1
            #elif KB == 66
               prefA(-128(pfA))
            #endif

            mulps 3024(pA), rB2
            addps rB2, rC22
            #if KB > 66
               movaps 3072(pB), rB2
            #elif KB == 66
               prefA((pfA))
            #endif

         #endif
         #if KB > 66
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3072(pA), rB0
            addps rB0, rC20
            #if KB > 67
               movaps 3088(pB), rB0
            #elif KB == 67
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3072(pA), rB1
            addps rB1, rC21
            #if KB > 67
               movaps 3104(pB), rB1
            #elif KB == 67
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 67
               movaps 3088(pA), rA0
            #elif KB == 67
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 67
               movaps 3104(pA), rA1
            #elif KB == 67
               prefA(-128(pfA))
            #endif

            mulps 3072(pA), rB2
            addps rB2, rC22
            #if KB > 67
               movaps 3120(pB), rB2
            #elif KB == 67
               prefA((pfA))
            #endif

         #endif
         #if KB > 67
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3120(pA), rB0
            addps rB0, rC20
            #if KB > 68
               movaps 3136(pB), rB0
            #elif KB == 68
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3120(pA), rB1
            addps rB1, rC21
            #if KB > 68
               movaps 3152(pB), rB1
            #elif KB == 68
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 68
               movaps 3136(pA), rA0
            #elif KB == 68
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 68
               movaps 3152(pA), rA1
            #elif KB == 68
               prefA(-128(pfA))
            #endif

            mulps 3120(pA), rB2
            addps rB2, rC22
            #if KB > 68
               movaps 3168(pB), rB2
            #elif KB == 68
               prefA((pfA))
            #endif

         #endif
         #if KB > 68
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3168(pA), rB0
            addps rB0, rC20
            #if KB > 69
               movaps 3184(pB), rB0
            #elif KB == 69
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3168(pA), rB1
            addps rB1, rC21
            #if KB > 69
               movaps 3200(pB), rB1
            #elif KB == 69
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 69
               movaps 3184(pA), rA0
            #elif KB == 69
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 69
               movaps 3200(pA), rA1
            #elif KB == 69
               prefA(-128(pfA))
            #endif

            mulps 3168(pA), rB2
            addps rB2, rC22
            #if KB > 69
               movaps 3216(pB), rB2
            #elif KB == 69
               prefA((pfA))
            #endif

         #endif
         #if KB > 69
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3216(pA), rB0
            addps rB0, rC20
            #if KB > 70
               movaps 3232(pB), rB0
            #elif KB == 70
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3216(pA), rB1
            addps rB1, rC21
            #if KB > 70
               movaps 3248(pB), rB1
            #elif KB == 70
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 70
               movaps 3232(pA), rA0
            #elif KB == 70
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 70
               movaps 3248(pA), rA1
            #elif KB == 70
               prefA(-128(pfA))
            #endif

            mulps 3216(pA), rB2
            addps rB2, rC22
            #if KB > 70
               movaps 3264(pB), rB2
            #elif KB == 70
               prefA((pfA))
            #endif

         #endif
         #if KB > 70
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3264(pA), rB0
            addps rB0, rC20
            #if KB > 71
               movaps 3280(pB), rB0
            #elif KB == 71
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3264(pA), rB1
            addps rB1, rC21
            #if KB > 71
               movaps 3296(pB), rB1
            #elif KB == 71
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 71
               movaps 3280(pA), rA0
            #elif KB == 71
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 71
               movaps 3296(pA), rA1
            #elif KB == 71
               prefA(-128(pfA))
            #endif

            mulps 3264(pA), rB2
            addps rB2, rC22
            #if KB > 71
               movaps 3312(pB), rB2
            #elif KB == 71
               prefA((pfA))
            #endif

         #endif
         #if KB > 71
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3312(pA), rB0
            addps rB0, rC20
            #if KB > 72
               movaps 3328(pB), rB0
            #elif KB == 72
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3312(pA), rB1
            addps rB1, rC21
            #if KB > 72
               movaps 3344(pB), rB1
            #elif KB == 72
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 72
               movaps 3328(pA), rA0
            #elif KB == 72
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 72
               movaps 3344(pA), rA1
            #elif KB == 72
               prefA(-128(pfA))
            #endif

            mulps 3312(pA), rB2
            addps rB2, rC22
            #if KB > 72
               movaps 3360(pB), rB2
            #elif KB == 72
               prefA((pfA))
            #endif

         #endif
         #if KB > 72
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3360(pA), rB0
            addps rB0, rC20
            #if KB > 73
               movaps 3376(pB), rB0
            #elif KB == 73
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3360(pA), rB1
            addps rB1, rC21
            #if KB > 73
               movaps 3392(pB), rB1
            #elif KB == 73
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 73
               movaps 3376(pA), rA0
            #elif KB == 73
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 73
               movaps 3392(pA), rA1
            #elif KB == 73
               prefA(-128(pfA))
            #endif

            mulps 3360(pA), rB2
            addps rB2, rC22
            #if KB > 73
               movaps 3408(pB), rB2
            #elif KB == 73
               prefA((pfA))
            #endif

         #endif
         #if KB > 73
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3408(pA), rB0
            addps rB0, rC20
            #if KB > 74
               movaps 3424(pB), rB0
            #elif KB == 74
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3408(pA), rB1
            addps rB1, rC21
            #if KB > 74
               movaps 3440(pB), rB1
            #elif KB == 74
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 74
               movaps 3424(pA), rA0
            #elif KB == 74
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 74
               movaps 3440(pA), rA1
            #elif KB == 74
               prefA(-128(pfA))
            #endif

            mulps 3408(pA), rB2
            addps rB2, rC22
            #if KB > 74
               movaps 3456(pB), rB2
            #elif KB == 74
               prefA((pfA))
            #endif

         #endif
         #if KB > 74
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3456(pA), rB0
            addps rB0, rC20
            #if KB > 75
               movaps 3472(pB), rB0
            #elif KB == 75
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3456(pA), rB1
            addps rB1, rC21
            #if KB > 75
               movaps 3488(pB), rB1
            #elif KB == 75
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 75
               movaps 3472(pA), rA0
            #elif KB == 75
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 75
               movaps 3488(pA), rA1
            #elif KB == 75
               prefA(-128(pfA))
            #endif

            mulps 3456(pA), rB2
            addps rB2, rC22
            #if KB > 75
               movaps 3504(pB), rB2
            #elif KB == 75
               prefA((pfA))
            #endif

         #endif
         #if KB > 75
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3504(pA), rB0
            addps rB0, rC20
            #if KB > 76
               movaps 3520(pB), rB0
            #elif KB == 76
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3504(pA), rB1
            addps rB1, rC21
            #if KB > 76
               movaps 3536(pB), rB1
            #elif KB == 76
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 76
               movaps 3520(pA), rA0
            #elif KB == 76
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 76
               movaps 3536(pA), rA1
            #elif KB == 76
               prefA(-128(pfA))
            #endif

            mulps 3504(pA), rB2
            addps rB2, rC22
            #if KB > 76
               movaps 3552(pB), rB2
            #elif KB == 76
               prefA((pfA))
            #endif

         #endif
         #if KB > 76
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3552(pA), rB0
            addps rB0, rC20
            #if KB > 77
               movaps 3568(pB), rB0
            #elif KB == 77
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3552(pA), rB1
            addps rB1, rC21
            #if KB > 77
               movaps 3584(pB), rB1
            #elif KB == 77
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 77
               movaps 3568(pA), rA0
            #elif KB == 77
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 77
               movaps 3584(pA), rA1
            #elif KB == 77
               prefA(-128(pfA))
            #endif

            mulps 3552(pA), rB2
            addps rB2, rC22
            #if KB > 77
               movaps 3600(pB), rB2
            #elif KB == 77
               prefA((pfA))
            #endif

         #endif
         #if KB > 77
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3600(pA), rB0
            addps rB0, rC20
            #if KB > 78
               movaps 3616(pB), rB0
            #elif KB == 78
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3600(pA), rB1
            addps rB1, rC21
            #if KB > 78
               movaps 3632(pB), rB1
            #elif KB == 78
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 78
               movaps 3616(pA), rA0
            #elif KB == 78
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 78
               movaps 3632(pA), rA1
            #elif KB == 78
               prefA(-128(pfA))
            #endif

            mulps 3600(pA), rB2
            addps rB2, rC22
            #if KB > 78
               movaps 3648(pB), rB2
            #elif KB == 78
               prefA((pfA))
            #endif

         #endif
         #if KB > 78
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3648(pA), rB0
            addps rB0, rC20
            #if KB > 79
               movaps 3664(pB), rB0
            #elif KB == 79
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3648(pA), rB1
            addps rB1, rC21
            #if KB > 79
               movaps 3680(pB), rB1
            #elif KB == 79
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 79
               movaps 3664(pA), rA0
            #elif KB == 79
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 79
               movaps 3680(pA), rA1
            #elif KB == 79
               prefA(-128(pfA))
            #endif

            mulps 3648(pA), rB2
            addps rB2, rC22
            #if KB > 79
               movaps 3696(pB), rB2
            #elif KB == 79
               prefA((pfA))
            #endif

         #endif
         #if KB > 79
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3696(pA), rB0
            addps rB0, rC20
            #if KB > 80
               movaps 3712(pB), rB0
            #elif KB == 80
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3696(pA), rB1
            addps rB1, rC21
            #if KB > 80
               movaps 3728(pB), rB1
            #elif KB == 80
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 80
               movaps 3712(pA), rA0
            #elif KB == 80
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 80
               movaps 3728(pA), rA1
            #elif KB == 80
               prefA(-128(pfA))
            #endif

            mulps 3696(pA), rB2
            addps rB2, rC22
            #if KB > 80
               movaps 3744(pB), rB2
            #elif KB == 80
               prefA((pfA))
            #endif

         #endif
         #if KB > 80
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3744(pA), rB0
            addps rB0, rC20
            #if KB > 81
               movaps 3760(pB), rB0
            #elif KB == 81
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3744(pA), rB1
            addps rB1, rC21
            #if KB > 81
               movaps 3776(pB), rB1
            #elif KB == 81
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 81
               movaps 3760(pA), rA0
            #elif KB == 81
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 81
               movaps 3776(pA), rA1
            #elif KB == 81
               prefA(-128(pfA))
            #endif

            mulps 3744(pA), rB2
            addps rB2, rC22
            #if KB > 81
               movaps 3792(pB), rB2
            #elif KB == 81
               prefA((pfA))
            #endif

         #endif
         #if KB > 81
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3792(pA), rB0
            addps rB0, rC20
            #if KB > 82
               movaps 3808(pB), rB0
            #elif KB == 82
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3792(pA), rB1
            addps rB1, rC21
            #if KB > 82
               movaps 3824(pB), rB1
            #elif KB == 82
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 82
               movaps 3808(pA), rA0
            #elif KB == 82
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 82
               movaps 3824(pA), rA1
            #elif KB == 82
               prefA(-128(pfA))
            #endif

            mulps 3792(pA), rB2
            addps rB2, rC22
            #if KB > 82
               movaps 3840(pB), rB2
            #elif KB == 82
               prefA((pfA))
            #endif

         #endif
         #if KB > 82
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3840(pA), rB0
            addps rB0, rC20
            #if KB > 83
               movaps 3856(pB), rB0
            #elif KB == 83
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3840(pA), rB1
            addps rB1, rC21
            #if KB > 83
               movaps 3872(pB), rB1
            #elif KB == 83
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 83
               movaps 3856(pA), rA0
            #elif KB == 83
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 83
               movaps 3872(pA), rA1
            #elif KB == 83
               prefA(-128(pfA))
            #endif

            mulps 3840(pA), rB2
            addps rB2, rC22
            #if KB > 83
               movaps 3888(pB), rB2
            #elif KB == 83
               prefA((pfA))
            #endif

         #endif
         #if KB > 83
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3888(pA), rB0
            addps rB0, rC20
            #if KB > 84
               movaps 3904(pB), rB0
            #elif KB == 84
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3888(pA), rB1
            addps rB1, rC21
            #if KB > 84
               movaps 3920(pB), rB1
            #elif KB == 84
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 84
               movaps 3904(pA), rA0
            #elif KB == 84
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 84
               movaps 3920(pA), rA1
            #elif KB == 84
               prefA(-128(pfA))
            #endif

            mulps 3888(pA), rB2
            addps rB2, rC22
            #if KB > 84
               movaps 3936(pB), rB2
            #elif KB == 84
               prefA((pfA))
            #endif

         #endif
         #if KB > 84
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3936(pA), rB0
            addps rB0, rC20
            #if KB > 85
               movaps 3952(pB), rB0
            #elif KB == 85
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3936(pA), rB1
            addps rB1, rC21
            #if KB > 85
               movaps 3968(pB), rB1
            #elif KB == 85
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 85
               movaps 3952(pA), rA0
            #elif KB == 85
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 85
               movaps 3968(pA), rA1
            #elif KB == 85
               prefA(-128(pfA))
            #endif

            mulps 3936(pA), rB2
            addps rB2, rC22
            #if KB > 85
               movaps 3984(pB), rB2
            #elif KB == 85
               prefA((pfA))
            #endif

         #endif
         #if KB > 85
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 3984(pA), rB0
            addps rB0, rC20
            #if KB > 86
               movaps 4000(pB), rB0
            #elif KB == 86
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 3984(pA), rB1
            addps rB1, rC21
            #if KB > 86
               movaps 4016(pB), rB1
            #elif KB == 86
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 86
               movaps 4000(pA), rA0
            #elif KB == 86
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 86
               movaps 4016(pA), rA1
            #elif KB == 86
               prefA(-128(pfA))
            #endif

            mulps 3984(pA), rB2
            addps rB2, rC22
            #if KB > 86
               movaps 4032(pB), rB2
            #elif KB == 86
               prefA((pfA))
            #endif

         #endif
         #if KB > 86
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4032(pA), rB0
            addps rB0, rC20
            #if KB > 87
               movaps 4048(pB), rB0
            #elif KB == 87
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4032(pA), rB1
            addps rB1, rC21
            #if KB > 87
               movaps 4064(pB), rB1
            #elif KB == 87
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 87
               movaps 4048(pA), rA0
            #elif KB == 87
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 87
               movaps 4064(pA), rA1
            #elif KB == 87
               prefA(-128(pfA))
            #endif

            mulps 4032(pA), rB2
            addps rB2, rC22
            #if KB > 87
               movaps 4080(pB), rB2
            #elif KB == 87
               prefA((pfA))
            #endif

         #endif
         #if KB > 87
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4080(pA), rB0
            addps rB0, rC20
            #if KB > 88
               movaps 4096(pB), rB0
            #elif KB == 88
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4080(pA), rB1
            addps rB1, rC21
            #if KB > 88
               movaps 4112(pB), rB1
            #elif KB == 88
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 88
               movaps 4096(pA), rA0
            #elif KB == 88
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 88
               movaps 4112(pA), rA1
            #elif KB == 88
               prefA(-128(pfA))
            #endif

            mulps 4080(pA), rB2
            addps rB2, rC22
            #if KB > 88
               movaps 4128(pB), rB2
            #elif KB == 88
               prefA((pfA))
            #endif

         #endif
         #if KB > 88
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4128(pA), rB0
            addps rB0, rC20
            #if KB > 89
               movaps 4144(pB), rB0
            #elif KB == 89
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4128(pA), rB1
            addps rB1, rC21
            #if KB > 89
               movaps 4160(pB), rB1
            #elif KB == 89
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 89
               movaps 4144(pA), rA0
            #elif KB == 89
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 89
               movaps 4160(pA), rA1
            #elif KB == 89
               prefA(-128(pfA))
            #endif

            mulps 4128(pA), rB2
            addps rB2, rC22
            #if KB > 89
               movaps 4176(pB), rB2
            #elif KB == 89
               prefA((pfA))
            #endif

         #endif
         #if KB > 89
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4176(pA), rB0
            addps rB0, rC20
            #if KB > 90
               movaps 4192(pB), rB0
            #elif KB == 90
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4176(pA), rB1
            addps rB1, rC21
            #if KB > 90
               movaps 4208(pB), rB1
            #elif KB == 90
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 90
               movaps 4192(pA), rA0
            #elif KB == 90
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 90
               movaps 4208(pA), rA1
            #elif KB == 90
               prefA(-128(pfA))
            #endif

            mulps 4176(pA), rB2
            addps rB2, rC22
            #if KB > 90
               movaps 4224(pB), rB2
            #elif KB == 90
               prefA((pfA))
            #endif

         #endif
         #if KB > 90
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4224(pA), rB0
            addps rB0, rC20
            #if KB > 91
               movaps 4240(pB), rB0
            #elif KB == 91
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4224(pA), rB1
            addps rB1, rC21
            #if KB > 91
               movaps 4256(pB), rB1
            #elif KB == 91
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 91
               movaps 4240(pA), rA0
            #elif KB == 91
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 91
               movaps 4256(pA), rA1
            #elif KB == 91
               prefA(-128(pfA))
            #endif

            mulps 4224(pA), rB2
            addps rB2, rC22
            #if KB > 91
               movaps 4272(pB), rB2
            #elif KB == 91
               prefA((pfA))
            #endif

         #endif
         #if KB > 91
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4272(pA), rB0
            addps rB0, rC20
            #if KB > 92
               movaps 4288(pB), rB0
            #elif KB == 92
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4272(pA), rB1
            addps rB1, rC21
            #if KB > 92
               movaps 4304(pB), rB1
            #elif KB == 92
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 92
               movaps 4288(pA), rA0
            #elif KB == 92
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 92
               movaps 4304(pA), rA1
            #elif KB == 92
               prefA(-128(pfA))
            #endif

            mulps 4272(pA), rB2
            addps rB2, rC22
            #if KB > 92
               movaps 4320(pB), rB2
            #elif KB == 92
               prefA((pfA))
            #endif

         #endif
         #if KB > 92
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4320(pA), rB0
            addps rB0, rC20
            #if KB > 93
               movaps 4336(pB), rB0
            #elif KB == 93
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4320(pA), rB1
            addps rB1, rC21
            #if KB > 93
               movaps 4352(pB), rB1
            #elif KB == 93
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 93
               movaps 4336(pA), rA0
            #elif KB == 93
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 93
               movaps 4352(pA), rA1
            #elif KB == 93
               prefA(-128(pfA))
            #endif

            mulps 4320(pA), rB2
            addps rB2, rC22
            #if KB > 93
               movaps 4368(pB), rB2
            #elif KB == 93
               prefA((pfA))
            #endif

         #endif
         #if KB > 93
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4368(pA), rB0
            addps rB0, rC20
            #if KB > 94
               movaps 4384(pB), rB0
            #elif KB == 94
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4368(pA), rB1
            addps rB1, rC21
            #if KB > 94
               movaps 4400(pB), rB1
            #elif KB == 94
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 94
               movaps 4384(pA), rA0
            #elif KB == 94
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 94
               movaps 4400(pA), rA1
            #elif KB == 94
               prefA(-128(pfA))
            #endif

            mulps 4368(pA), rB2
            addps rB2, rC22
            #if KB > 94
               movaps 4416(pB), rB2
            #elif KB == 94
               prefA((pfA))
            #endif

         #endif
         #if KB > 94
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4416(pA), rB0
            addps rB0, rC20
            #if KB > 95
               movaps 4432(pB), rB0
            #elif KB == 95
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4416(pA), rB1
            addps rB1, rC21
            #if KB > 95
               movaps 4448(pB), rB1
            #elif KB == 95
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 95
               movaps 4432(pA), rA0
            #elif KB == 95
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 95
               movaps 4448(pA), rA1
            #elif KB == 95
               prefA(-128(pfA))
            #endif

            mulps 4416(pA), rB2
            addps rB2, rC22
            #if KB > 95
               movaps 4464(pB), rB2
            #elif KB == 95
               prefA((pfA))
            #endif

         #endif
         #if KB > 95
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4464(pA), rB0
            addps rB0, rC20
            #if KB > 96
               movaps 4480(pB), rB0
            #elif KB == 96
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4464(pA), rB1
            addps rB1, rC21
            #if KB > 96
               movaps 4496(pB), rB1
            #elif KB == 96
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 96
               movaps 4480(pA), rA0
            #elif KB == 96
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 96
               movaps 4496(pA), rA1
            #elif KB == 96
               prefA(-128(pfA))
            #endif

            mulps 4464(pA), rB2
            addps rB2, rC22
            #if KB > 96
               movaps 4512(pB), rB2
            #elif KB == 96
               prefA((pfA))
            #endif

         #endif
         #if KB > 96
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4512(pA), rB0
            addps rB0, rC20
            #if KB > 97
               movaps 4528(pB), rB0
            #elif KB == 97
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4512(pA), rB1
            addps rB1, rC21
            #if KB > 97
               movaps 4544(pB), rB1
            #elif KB == 97
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 97
               movaps 4528(pA), rA0
            #elif KB == 97
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 97
               movaps 4544(pA), rA1
            #elif KB == 97
               prefA(-128(pfA))
            #endif

            mulps 4512(pA), rB2
            addps rB2, rC22
            #if KB > 97
               movaps 4560(pB), rB2
            #elif KB == 97
               prefA((pfA))
            #endif

         #endif
         #if KB > 97
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4560(pA), rB0
            addps rB0, rC20
            #if KB > 98
               movaps 4576(pB), rB0
            #elif KB == 98
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4560(pA), rB1
            addps rB1, rC21
            #if KB > 98
               movaps 4592(pB), rB1
            #elif KB == 98
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 98
               movaps 4576(pA), rA0
            #elif KB == 98
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 98
               movaps 4592(pA), rA1
            #elif KB == 98
               prefA(-128(pfA))
            #endif

            mulps 4560(pA), rB2
            addps rB2, rC22
            #if KB > 98
               movaps 4608(pB), rB2
            #elif KB == 98
               prefA((pfA))
            #endif

         #endif
         #if KB > 98
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4608(pA), rB0
            addps rB0, rC20
            #if KB > 99
               movaps 4624(pB), rB0
            #elif KB == 99
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4608(pA), rB1
            addps rB1, rC21
            #if KB > 99
               movaps 4640(pB), rB1
            #elif KB == 99
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 99
               movaps 4624(pA), rA0
            #elif KB == 99
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 99
               movaps 4640(pA), rA1
            #elif KB == 99
               prefA(-128(pfA))
            #endif

            mulps 4608(pA), rB2
            addps rB2, rC22
            #if KB > 99
               movaps 4656(pB), rB2
            #elif KB == 99
               prefA((pfA))
            #endif

         #endif
         #if KB > 99
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4656(pA), rB0
            addps rB0, rC20
            #if KB > 100
               movaps 4672(pB), rB0
            #elif KB == 100
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4656(pA), rB1
            addps rB1, rC21
            #if KB > 100
               movaps 4688(pB), rB1
            #elif KB == 100
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 100
               movaps 4672(pA), rA0
            #elif KB == 100
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 100
               movaps 4688(pA), rA1
            #elif KB == 100
               prefA(-128(pfA))
            #endif

            mulps 4656(pA), rB2
            addps rB2, rC22
            #if KB > 100
               movaps 4704(pB), rB2
            #elif KB == 100
               prefA((pfA))
            #endif

         #endif
         #if KB > 100
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4704(pA), rB0
            addps rB0, rC20
            #if KB > 101
               movaps 4720(pB), rB0
            #elif KB == 101
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4704(pA), rB1
            addps rB1, rC21
            #if KB > 101
               movaps 4736(pB), rB1
            #elif KB == 101
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 101
               movaps 4720(pA), rA0
            #elif KB == 101
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 101
               movaps 4736(pA), rA1
            #elif KB == 101
               prefA(-128(pfA))
            #endif

            mulps 4704(pA), rB2
            addps rB2, rC22
            #if KB > 101
               movaps 4752(pB), rB2
            #elif KB == 101
               prefA((pfA))
            #endif

         #endif
         #if KB > 101
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4752(pA), rB0
            addps rB0, rC20
            #if KB > 102
               movaps 4768(pB), rB0
            #elif KB == 102
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4752(pA), rB1
            addps rB1, rC21
            #if KB > 102
               movaps 4784(pB), rB1
            #elif KB == 102
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 102
               movaps 4768(pA), rA0
            #elif KB == 102
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 102
               movaps 4784(pA), rA1
            #elif KB == 102
               prefA(-128(pfA))
            #endif

            mulps 4752(pA), rB2
            addps rB2, rC22
            #if KB > 102
               movaps 4800(pB), rB2
            #elif KB == 102
               prefA((pfA))
            #endif

         #endif
         #if KB > 102
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4800(pA), rB0
            addps rB0, rC20
            #if KB > 103
               movaps 4816(pB), rB0
            #elif KB == 103
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4800(pA), rB1
            addps rB1, rC21
            #if KB > 103
               movaps 4832(pB), rB1
            #elif KB == 103
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 103
               movaps 4816(pA), rA0
            #elif KB == 103
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 103
               movaps 4832(pA), rA1
            #elif KB == 103
               prefA(-128(pfA))
            #endif

            mulps 4800(pA), rB2
            addps rB2, rC22
            #if KB > 103
               movaps 4848(pB), rB2
            #elif KB == 103
               prefA((pfA))
            #endif

         #endif
         #if KB > 103
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4848(pA), rB0
            addps rB0, rC20
            #if KB > 104
               movaps 4864(pB), rB0
            #elif KB == 104
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4848(pA), rB1
            addps rB1, rC21
            #if KB > 104
               movaps 4880(pB), rB1
            #elif KB == 104
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 104
               movaps 4864(pA), rA0
            #elif KB == 104
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 104
               movaps 4880(pA), rA1
            #elif KB == 104
               prefA(-128(pfA))
            #endif

            mulps 4848(pA), rB2
            addps rB2, rC22
            #if KB > 104
               movaps 4896(pB), rB2
            #elif KB == 104
               prefA((pfA))
            #endif

         #endif
         #if KB > 104
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4896(pA), rB0
            addps rB0, rC20
            #if KB > 105
               movaps 4912(pB), rB0
            #elif KB == 105
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4896(pA), rB1
            addps rB1, rC21
            #if KB > 105
               movaps 4928(pB), rB1
            #elif KB == 105
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 105
               movaps 4912(pA), rA0
            #elif KB == 105
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 105
               movaps 4928(pA), rA1
            #elif KB == 105
               prefA(-128(pfA))
            #endif

            mulps 4896(pA), rB2
            addps rB2, rC22
            #if KB > 105
               movaps 4944(pB), rB2
            #elif KB == 105
               prefA((pfA))
            #endif

         #endif
         #if KB > 105
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4944(pA), rB0
            addps rB0, rC20
            #if KB > 106
               movaps 4960(pB), rB0
            #elif KB == 106
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4944(pA), rB1
            addps rB1, rC21
            #if KB > 106
               movaps 4976(pB), rB1
            #elif KB == 106
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 106
               movaps 4960(pA), rA0
            #elif KB == 106
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 106
               movaps 4976(pA), rA1
            #elif KB == 106
               prefA(-128(pfA))
            #endif

            mulps 4944(pA), rB2
            addps rB2, rC22
            #if KB > 106
               movaps 4992(pB), rB2
            #elif KB == 106
               prefA((pfA))
            #endif

         #endif
         #if KB > 106
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 4992(pA), rB0
            addps rB0, rC20
            #if KB > 107
               movaps 5008(pB), rB0
            #elif KB == 107
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 4992(pA), rB1
            addps rB1, rC21
            #if KB > 107
               movaps 5024(pB), rB1
            #elif KB == 107
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 107
               movaps 5008(pA), rA0
            #elif KB == 107
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 107
               movaps 5024(pA), rA1
            #elif KB == 107
               prefA(-128(pfA))
            #endif

            mulps 4992(pA), rB2
            addps rB2, rC22
            #if KB > 107
               movaps 5040(pB), rB2
            #elif KB == 107
               prefA((pfA))
            #endif

         #endif
         #if KB > 107
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5040(pA), rB0
            addps rB0, rC20
            #if KB > 108
               movaps 5056(pB), rB0
            #elif KB == 108
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5040(pA), rB1
            addps rB1, rC21
            #if KB > 108
               movaps 5072(pB), rB1
            #elif KB == 108
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 108
               movaps 5056(pA), rA0
            #elif KB == 108
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 108
               movaps 5072(pA), rA1
            #elif KB == 108
               prefA(-128(pfA))
            #endif

            mulps 5040(pA), rB2
            addps rB2, rC22
            #if KB > 108
               movaps 5088(pB), rB2
            #elif KB == 108
               prefA((pfA))
            #endif

         #endif
         #if KB > 108
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5088(pA), rB0
            addps rB0, rC20
            #if KB > 109
               movaps 5104(pB), rB0
            #elif KB == 109
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5088(pA), rB1
            addps rB1, rC21
            #if KB > 109
               movaps 5120(pB), rB1
            #elif KB == 109
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 109
               movaps 5104(pA), rA0
            #elif KB == 109
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 109
               movaps 5120(pA), rA1
            #elif KB == 109
               prefA(-128(pfA))
            #endif

            mulps 5088(pA), rB2
            addps rB2, rC22
            #if KB > 109
               movaps 5136(pB), rB2
            #elif KB == 109
               prefA((pfA))
            #endif

         #endif
         #if KB > 109
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5136(pA), rB0
            addps rB0, rC20
            #if KB > 110
               movaps 5152(pB), rB0
            #elif KB == 110
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5136(pA), rB1
            addps rB1, rC21
            #if KB > 110
               movaps 5168(pB), rB1
            #elif KB == 110
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 110
               movaps 5152(pA), rA0
            #elif KB == 110
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 110
               movaps 5168(pA), rA1
            #elif KB == 110
               prefA(-128(pfA))
            #endif

            mulps 5136(pA), rB2
            addps rB2, rC22
            #if KB > 110
               movaps 5184(pB), rB2
            #elif KB == 110
               prefA((pfA))
            #endif

         #endif
         #if KB > 110
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5184(pA), rB0
            addps rB0, rC20
            #if KB > 111
               movaps 5200(pB), rB0
            #elif KB == 111
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5184(pA), rB1
            addps rB1, rC21
            #if KB > 111
               movaps 5216(pB), rB1
            #elif KB == 111
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 111
               movaps 5200(pA), rA0
            #elif KB == 111
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 111
               movaps 5216(pA), rA1
            #elif KB == 111
               prefA(-128(pfA))
            #endif

            mulps 5184(pA), rB2
            addps rB2, rC22
            #if KB > 111
               movaps 5232(pB), rB2
            #elif KB == 111
               prefA((pfA))
            #endif

         #endif
         #if KB > 111
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5232(pA), rB0
            addps rB0, rC20
            #if KB > 112
               movaps 5248(pB), rB0
            #elif KB == 112
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5232(pA), rB1
            addps rB1, rC21
            #if KB > 112
               movaps 5264(pB), rB1
            #elif KB == 112
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 112
               movaps 5248(pA), rA0
            #elif KB == 112
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 112
               movaps 5264(pA), rA1
            #elif KB == 112
               prefA(-128(pfA))
            #endif

            mulps 5232(pA), rB2
            addps rB2, rC22
            #if KB > 112
               movaps 5280(pB), rB2
            #elif KB == 112
               prefA((pfA))
            #endif

         #endif
         #if KB > 112
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5280(pA), rB0
            addps rB0, rC20
            #if KB > 113
               movaps 5296(pB), rB0
            #elif KB == 113
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5280(pA), rB1
            addps rB1, rC21
            #if KB > 113
               movaps 5312(pB), rB1
            #elif KB == 113
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 113
               movaps 5296(pA), rA0
            #elif KB == 113
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 113
               movaps 5312(pA), rA1
            #elif KB == 113
               prefA(-128(pfA))
            #endif

            mulps 5280(pA), rB2
            addps rB2, rC22
            #if KB > 113
               movaps 5328(pB), rB2
            #elif KB == 113
               prefA((pfA))
            #endif

         #endif
         #if KB > 113
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5328(pA), rB0
            addps rB0, rC20
            #if KB > 114
               movaps 5344(pB), rB0
            #elif KB == 114
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5328(pA), rB1
            addps rB1, rC21
            #if KB > 114
               movaps 5360(pB), rB1
            #elif KB == 114
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 114
               movaps 5344(pA), rA0
            #elif KB == 114
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 114
               movaps 5360(pA), rA1
            #elif KB == 114
               prefA(-128(pfA))
            #endif

            mulps 5328(pA), rB2
            addps rB2, rC22
            #if KB > 114
               movaps 5376(pB), rB2
            #elif KB == 114
               prefA((pfA))
            #endif

         #endif
         #if KB > 114
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5376(pA), rB0
            addps rB0, rC20
            #if KB > 115
               movaps 5392(pB), rB0
            #elif KB == 115
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5376(pA), rB1
            addps rB1, rC21
            #if KB > 115
               movaps 5408(pB), rB1
            #elif KB == 115
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 115
               movaps 5392(pA), rA0
            #elif KB == 115
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 115
               movaps 5408(pA), rA1
            #elif KB == 115
               prefA(-128(pfA))
            #endif

            mulps 5376(pA), rB2
            addps rB2, rC22
            #if KB > 115
               movaps 5424(pB), rB2
            #elif KB == 115
               prefA((pfA))
            #endif

         #endif
         #if KB > 115
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5424(pA), rB0
            addps rB0, rC20
            #if KB > 116
               movaps 5440(pB), rB0
            #elif KB == 116
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5424(pA), rB1
            addps rB1, rC21
            #if KB > 116
               movaps 5456(pB), rB1
            #elif KB == 116
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 116
               movaps 5440(pA), rA0
            #elif KB == 116
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 116
               movaps 5456(pA), rA1
            #elif KB == 116
               prefA(-128(pfA))
            #endif

            mulps 5424(pA), rB2
            addps rB2, rC22
            #if KB > 116
               movaps 5472(pB), rB2
            #elif KB == 116
               prefA((pfA))
            #endif

         #endif
         #if KB > 116
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5472(pA), rB0
            addps rB0, rC20
            #if KB > 117
               movaps 5488(pB), rB0
            #elif KB == 117
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5472(pA), rB1
            addps rB1, rC21
            #if KB > 117
               movaps 5504(pB), rB1
            #elif KB == 117
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 117
               movaps 5488(pA), rA0
            #elif KB == 117
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 117
               movaps 5504(pA), rA1
            #elif KB == 117
               prefA(-128(pfA))
            #endif

            mulps 5472(pA), rB2
            addps rB2, rC22
            #if KB > 117
               movaps 5520(pB), rB2
            #elif KB == 117
               prefA((pfA))
            #endif

         #endif
         #if KB > 117
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5520(pA), rB0
            addps rB0, rC20
            #if KB > 118
               movaps 5536(pB), rB0
            #elif KB == 118
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5520(pA), rB1
            addps rB1, rC21
            #if KB > 118
               movaps 5552(pB), rB1
            #elif KB == 118
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 118
               movaps 5536(pA), rA0
            #elif KB == 118
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 118
               movaps 5552(pA), rA1
            #elif KB == 118
               prefA(-128(pfA))
            #endif

            mulps 5520(pA), rB2
            addps rB2, rC22
            #if KB > 118
               movaps 5568(pB), rB2
            #elif KB == 118
               prefA((pfA))
            #endif

         #endif
         #if KB > 118
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5568(pA), rB0
            addps rB0, rC20
            #if KB > 119
               movaps 5584(pB), rB0
            #elif KB == 119
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5568(pA), rB1
            addps rB1, rC21
            #if KB > 119
               movaps 5600(pB), rB1
            #elif KB == 119
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 119
               movaps 5584(pA), rA0
            #elif KB == 119
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 119
               movaps 5600(pA), rA1
            #elif KB == 119
               prefA(-128(pfA))
            #endif

            mulps 5568(pA), rB2
            addps rB2, rC22
            #if KB > 119
               movaps 5616(pB), rB2
            #elif KB == 119
               prefA((pfA))
            #endif

         #endif
         #if KB > 119
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5616(pA), rB0
            addps rB0, rC20
            #if KB > 120
               movaps 5632(pB), rB0
            #elif KB == 120
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5616(pA), rB1
            addps rB1, rC21
            #if KB > 120
               movaps 5648(pB), rB1
            #elif KB == 120
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 120
               movaps 5632(pA), rA0
            #elif KB == 120
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 120
               movaps 5648(pA), rA1
            #elif KB == 120
               prefA(-128(pfA))
            #endif

            mulps 5616(pA), rB2
            addps rB2, rC22
            #if KB > 120
               movaps 5664(pB), rB2
            #elif KB == 120
               prefA((pfA))
            #endif

         #endif
         #if KB > 120
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5664(pA), rB0
            addps rB0, rC20
            #if KB > 121
               movaps 5680(pB), rB0
            #elif KB == 121
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5664(pA), rB1
            addps rB1, rC21
            #if KB > 121
               movaps 5696(pB), rB1
            #elif KB == 121
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 121
               movaps 5680(pA), rA0
            #elif KB == 121
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 121
               movaps 5696(pA), rA1
            #elif KB == 121
               prefA(-128(pfA))
            #endif

            mulps 5664(pA), rB2
            addps rB2, rC22
            #if KB > 121
               movaps 5712(pB), rB2
            #elif KB == 121
               prefA((pfA))
            #endif

         #endif
         #if KB > 121
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5712(pA), rB0
            addps rB0, rC20
            #if KB > 122
               movaps 5728(pB), rB0
            #elif KB == 122
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5712(pA), rB1
            addps rB1, rC21
            #if KB > 122
               movaps 5744(pB), rB1
            #elif KB == 122
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 122
               movaps 5728(pA), rA0
            #elif KB == 122
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 122
               movaps 5744(pA), rA1
            #elif KB == 122
               prefA(-128(pfA))
            #endif

            mulps 5712(pA), rB2
            addps rB2, rC22
            #if KB > 122
               movaps 5760(pB), rB2
            #elif KB == 122
               prefA((pfA))
            #endif

         #endif
         #if KB > 122
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5760(pA), rB0
            addps rB0, rC20
            #if KB > 123
               movaps 5776(pB), rB0
            #elif KB == 123
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5760(pA), rB1
            addps rB1, rC21
            #if KB > 123
               movaps 5792(pB), rB1
            #elif KB == 123
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 123
               movaps 5776(pA), rA0
            #elif KB == 123
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 123
               movaps 5792(pA), rA1
            #elif KB == 123
               prefA(-128(pfA))
            #endif

            mulps 5760(pA), rB2
            addps rB2, rC22
            #if KB > 123
               movaps 5808(pB), rB2
            #elif KB == 123
               prefA((pfA))
            #endif

         #endif
         #if KB > 123
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5808(pA), rB0
            addps rB0, rC20
            #if KB > 124
               movaps 5824(pB), rB0
            #elif KB == 124
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5808(pA), rB1
            addps rB1, rC21
            #if KB > 124
               movaps 5840(pB), rB1
            #elif KB == 124
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 124
               movaps 5824(pA), rA0
            #elif KB == 124
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 124
               movaps 5840(pA), rA1
            #elif KB == 124
               prefA(-128(pfA))
            #endif

            mulps 5808(pA), rB2
            addps rB2, rC22
            #if KB > 124
               movaps 5856(pB), rB2
            #elif KB == 124
               prefA((pfA))
            #endif

         #endif
         #if KB > 124
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5856(pA), rB0
            addps rB0, rC20
            #if KB > 125
               movaps 5872(pB), rB0
            #elif KB == 125
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5856(pA), rB1
            addps rB1, rC21
            #if KB > 125
               movaps 5888(pB), rB1
            #elif KB == 125
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 125
               movaps 5872(pA), rA0
            #elif KB == 125
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 125
               movaps 5888(pA), rA1
            #elif KB == 125
               prefA(-128(pfA))
            #endif

            mulps 5856(pA), rB2
            addps rB2, rC22
            #if KB > 125
               movaps 5904(pB), rB2
            #elif KB == 125
               prefA((pfA))
            #endif

         #endif
         #if KB > 125
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5904(pA), rB0
            addps rB0, rC20
            #if KB > 126
               movaps 5920(pB), rB0
            #elif KB == 126
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5904(pA), rB1
            addps rB1, rC21
            #if KB > 126
               movaps 5936(pB), rB1
            #elif KB == 126
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 126
               movaps 5920(pA), rA0
            #elif KB == 126
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 126
               movaps 5936(pA), rA1
            #elif KB == 126
               prefA(-128(pfA))
            #endif

            mulps 5904(pA), rB2
            addps rB2, rC22
            #if KB > 126
               movaps 5952(pB), rB2
            #elif KB == 126
               prefA((pfA))
            #endif

         #endif
         #if KB > 126
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 5952(pA), rB0
            addps rB0, rC20
            #if KB > 127
               movaps 5968(pB), rB0
            #elif KB == 127
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 5952(pA), rB1
            addps rB1, rC21
            #if KB > 127
               movaps 5984(pB), rB1
            #elif KB == 127
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 127
               movaps 5968(pA), rA0
            #elif KB == 127
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 127
               movaps 5984(pA), rA1
            #elif KB == 127
               prefA(-128(pfA))
            #endif

            mulps 5952(pA), rB2
            addps rB2, rC22
            #if KB > 127
               movaps 6000(pB), rB2
            #elif KB == 127
               prefA((pfA))
            #endif

         #endif
         #if KB > 127
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6000(pA), rB0
            addps rB0, rC20
            #if KB > 128
               movaps 6016(pB), rB0
            #elif KB == 128
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6000(pA), rB1
            addps rB1, rC21
            #if KB > 128
               movaps 6032(pB), rB1
            #elif KB == 128
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 128
               movaps 6016(pA), rA0
            #elif KB == 128
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 128
               movaps 6032(pA), rA1
            #elif KB == 128
               prefA(-128(pfA))
            #endif

            mulps 6000(pA), rB2
            addps rB2, rC22
            #if KB > 128
               movaps 6048(pB), rB2
            #elif KB == 128
               prefA((pfA))
            #endif

         #endif
         #if KB > 128
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6048(pA), rB0
            addps rB0, rC20
            #if KB > 129
               movaps 6064(pB), rB0
            #elif KB == 129
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6048(pA), rB1
            addps rB1, rC21
            #if KB > 129
               movaps 6080(pB), rB1
            #elif KB == 129
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 129
               movaps 6064(pA), rA0
            #elif KB == 129
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 129
               movaps 6080(pA), rA1
            #elif KB == 129
               prefA(-128(pfA))
            #endif

            mulps 6048(pA), rB2
            addps rB2, rC22
            #if KB > 129
               movaps 6096(pB), rB2
            #elif KB == 129
               prefA((pfA))
            #endif

         #endif
         #if KB > 129
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6096(pA), rB0
            addps rB0, rC20
            #if KB > 130
               movaps 6112(pB), rB0
            #elif KB == 130
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6096(pA), rB1
            addps rB1, rC21
            #if KB > 130
               movaps 6128(pB), rB1
            #elif KB == 130
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 130
               movaps 6112(pA), rA0
            #elif KB == 130
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 130
               movaps 6128(pA), rA1
            #elif KB == 130
               prefA(-128(pfA))
            #endif

            mulps 6096(pA), rB2
            addps rB2, rC22
            #if KB > 130
               movaps 6144(pB), rB2
            #elif KB == 130
               prefA((pfA))
            #endif

         #endif
         #if KB > 130
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6144(pA), rB0
            addps rB0, rC20
            #if KB > 131
               movaps 6160(pB), rB0
            #elif KB == 131
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6144(pA), rB1
            addps rB1, rC21
            #if KB > 131
               movaps 6176(pB), rB1
            #elif KB == 131
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 131
               movaps 6160(pA), rA0
            #elif KB == 131
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 131
               movaps 6176(pA), rA1
            #elif KB == 131
               prefA(-128(pfA))
            #endif

            mulps 6144(pA), rB2
            addps rB2, rC22
            #if KB > 131
               movaps 6192(pB), rB2
            #elif KB == 131
               prefA((pfA))
            #endif

         #endif
         #if KB > 131
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6192(pA), rB0
            addps rB0, rC20
            #if KB > 132
               movaps 6208(pB), rB0
            #elif KB == 132
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6192(pA), rB1
            addps rB1, rC21
            #if KB > 132
               movaps 6224(pB), rB1
            #elif KB == 132
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 132
               movaps 6208(pA), rA0
            #elif KB == 132
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 132
               movaps 6224(pA), rA1
            #elif KB == 132
               prefA(-128(pfA))
            #endif

            mulps 6192(pA), rB2
            addps rB2, rC22
            #if KB > 132
               movaps 6240(pB), rB2
            #elif KB == 132
               prefA((pfA))
            #endif

         #endif
         #if KB > 132
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6240(pA), rB0
            addps rB0, rC20
            #if KB > 133
               movaps 6256(pB), rB0
            #elif KB == 133
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6240(pA), rB1
            addps rB1, rC21
            #if KB > 133
               movaps 6272(pB), rB1
            #elif KB == 133
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 133
               movaps 6256(pA), rA0
            #elif KB == 133
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 133
               movaps 6272(pA), rA1
            #elif KB == 133
               prefA(-128(pfA))
            #endif

            mulps 6240(pA), rB2
            addps rB2, rC22
            #if KB > 133
               movaps 6288(pB), rB2
            #elif KB == 133
               prefA((pfA))
            #endif

         #endif
         #if KB > 133
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6288(pA), rB0
            addps rB0, rC20
            #if KB > 134
               movaps 6304(pB), rB0
            #elif KB == 134
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6288(pA), rB1
            addps rB1, rC21
            #if KB > 134
               movaps 6320(pB), rB1
            #elif KB == 134
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 134
               movaps 6304(pA), rA0
            #elif KB == 134
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 134
               movaps 6320(pA), rA1
            #elif KB == 134
               prefA(-128(pfA))
            #endif

            mulps 6288(pA), rB2
            addps rB2, rC22
            #if KB > 134
               movaps 6336(pB), rB2
            #elif KB == 134
               prefA((pfA))
            #endif

         #endif
         #if KB > 134
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6336(pA), rB0
            addps rB0, rC20
            #if KB > 135
               movaps 6352(pB), rB0
            #elif KB == 135
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6336(pA), rB1
            addps rB1, rC21
            #if KB > 135
               movaps 6368(pB), rB1
            #elif KB == 135
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 135
               movaps 6352(pA), rA0
            #elif KB == 135
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 135
               movaps 6368(pA), rA1
            #elif KB == 135
               prefA(-128(pfA))
            #endif

            mulps 6336(pA), rB2
            addps rB2, rC22
            #if KB > 135
               movaps 6384(pB), rB2
            #elif KB == 135
               prefA((pfA))
            #endif

         #endif
         #if KB > 135
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6384(pA), rB0
            addps rB0, rC20
            #if KB > 136
               movaps 6400(pB), rB0
            #elif KB == 136
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6384(pA), rB1
            addps rB1, rC21
            #if KB > 136
               movaps 6416(pB), rB1
            #elif KB == 136
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 136
               movaps 6400(pA), rA0
            #elif KB == 136
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 136
               movaps 6416(pA), rA1
            #elif KB == 136
               prefA(-128(pfA))
            #endif

            mulps 6384(pA), rB2
            addps rB2, rC22
            #if KB > 136
               movaps 6432(pB), rB2
            #elif KB == 136
               prefA((pfA))
            #endif

         #endif
         #if KB > 136
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6432(pA), rB0
            addps rB0, rC20
            #if KB > 137
               movaps 6448(pB), rB0
            #elif KB == 137
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6432(pA), rB1
            addps rB1, rC21
            #if KB > 137
               movaps 6464(pB), rB1
            #elif KB == 137
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 137
               movaps 6448(pA), rA0
            #elif KB == 137
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 137
               movaps 6464(pA), rA1
            #elif KB == 137
               prefA(-128(pfA))
            #endif

            mulps 6432(pA), rB2
            addps rB2, rC22
            #if KB > 137
               movaps 6480(pB), rB2
            #elif KB == 137
               prefA((pfA))
            #endif

         #endif
         #if KB > 137
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6480(pA), rB0
            addps rB0, rC20
            #if KB > 138
               movaps 6496(pB), rB0
            #elif KB == 138
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6480(pA), rB1
            addps rB1, rC21
            #if KB > 138
               movaps 6512(pB), rB1
            #elif KB == 138
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 138
               movaps 6496(pA), rA0
            #elif KB == 138
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 138
               movaps 6512(pA), rA1
            #elif KB == 138
               prefA(-128(pfA))
            #endif

            mulps 6480(pA), rB2
            addps rB2, rC22
            #if KB > 138
               movaps 6528(pB), rB2
            #elif KB == 138
               prefA((pfA))
            #endif

         #endif
         #if KB > 138
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6528(pA), rB0
            addps rB0, rC20
            #if KB > 139
               movaps 6544(pB), rB0
            #elif KB == 139
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6528(pA), rB1
            addps rB1, rC21
            #if KB > 139
               movaps 6560(pB), rB1
            #elif KB == 139
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 139
               movaps 6544(pA), rA0
            #elif KB == 139
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 139
               movaps 6560(pA), rA1
            #elif KB == 139
               prefA(-128(pfA))
            #endif

            mulps 6528(pA), rB2
            addps rB2, rC22
            #if KB > 139
               movaps 6576(pB), rB2
            #elif KB == 139
               prefA((pfA))
            #endif

         #endif
         #if KB > 139
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6576(pA), rB0
            addps rB0, rC20
            #if KB > 140
               movaps 6592(pB), rB0
            #elif KB == 140
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6576(pA), rB1
            addps rB1, rC21
            #if KB > 140
               movaps 6608(pB), rB1
            #elif KB == 140
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 140
               movaps 6592(pA), rA0
            #elif KB == 140
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 140
               movaps 6608(pA), rA1
            #elif KB == 140
               prefA(-128(pfA))
            #endif

            mulps 6576(pA), rB2
            addps rB2, rC22
            #if KB > 140
               movaps 6624(pB), rB2
            #elif KB == 140
               prefA((pfA))
            #endif

         #endif
         #if KB > 140
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6624(pA), rB0
            addps rB0, rC20
            #if KB > 141
               movaps 6640(pB), rB0
            #elif KB == 141
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6624(pA), rB1
            addps rB1, rC21
            #if KB > 141
               movaps 6656(pB), rB1
            #elif KB == 141
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 141
               movaps 6640(pA), rA0
            #elif KB == 141
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 141
               movaps 6656(pA), rA1
            #elif KB == 141
               prefA(-128(pfA))
            #endif

            mulps 6624(pA), rB2
            addps rB2, rC22
            #if KB > 141
               movaps 6672(pB), rB2
            #elif KB == 141
               prefA((pfA))
            #endif

         #endif
         #if KB > 141
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6672(pA), rB0
            addps rB0, rC20
            #if KB > 142
               movaps 6688(pB), rB0
            #elif KB == 142
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6672(pA), rB1
            addps rB1, rC21
            #if KB > 142
               movaps 6704(pB), rB1
            #elif KB == 142
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 142
               movaps 6688(pA), rA0
            #elif KB == 142
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 142
               movaps 6704(pA), rA1
            #elif KB == 142
               prefA(-128(pfA))
            #endif

            mulps 6672(pA), rB2
            addps rB2, rC22
            #if KB > 142
               movaps 6720(pB), rB2
            #elif KB == 142
               prefA((pfA))
            #endif

         #endif
         #if KB > 142
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6720(pA), rB0
            addps rB0, rC20
            #if KB > 143
               movaps 6736(pB), rB0
            #elif KB == 143
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6720(pA), rB1
            addps rB1, rC21
            #if KB > 143
               movaps 6752(pB), rB1
            #elif KB == 143
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 143
               movaps 6736(pA), rA0
            #elif KB == 143
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 143
               movaps 6752(pA), rA1
            #elif KB == 143
               prefA(-128(pfA))
            #endif

            mulps 6720(pA), rB2
            addps rB2, rC22
            #if KB > 143
               movaps 6768(pB), rB2
            #elif KB == 143
               prefA((pfA))
            #endif

         #endif
         #if KB > 143
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6768(pA), rB0
            addps rB0, rC20
            #if KB > 144
               movaps 6784(pB), rB0
            #elif KB == 144
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6768(pA), rB1
            addps rB1, rC21
            #if KB > 144
               movaps 6800(pB), rB1
            #elif KB == 144
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 144
               movaps 6784(pA), rA0
            #elif KB == 144
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 144
               movaps 6800(pA), rA1
            #elif KB == 144
               prefA(-128(pfA))
            #endif

            mulps 6768(pA), rB2
            addps rB2, rC22
            #if KB > 144
               movaps 6816(pB), rB2
            #elif KB == 144
               prefA((pfA))
            #endif

         #endif
         #if KB > 144
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6816(pA), rB0
            addps rB0, rC20
            #if KB > 145
               movaps 6832(pB), rB0
            #elif KB == 145
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6816(pA), rB1
            addps rB1, rC21
            #if KB > 145
               movaps 6848(pB), rB1
            #elif KB == 145
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 145
               movaps 6832(pA), rA0
            #elif KB == 145
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 145
               movaps 6848(pA), rA1
            #elif KB == 145
               prefA(-128(pfA))
            #endif

            mulps 6816(pA), rB2
            addps rB2, rC22
            #if KB > 145
               movaps 6864(pB), rB2
            #elif KB == 145
               prefA((pfA))
            #endif

         #endif
         #if KB > 145
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6864(pA), rB0
            addps rB0, rC20
            #if KB > 146
               movaps 6880(pB), rB0
            #elif KB == 146
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6864(pA), rB1
            addps rB1, rC21
            #if KB > 146
               movaps 6896(pB), rB1
            #elif KB == 146
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 146
               movaps 6880(pA), rA0
            #elif KB == 146
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 146
               movaps 6896(pA), rA1
            #elif KB == 146
               prefA(-128(pfA))
            #endif

            mulps 6864(pA), rB2
            addps rB2, rC22
            #if KB > 146
               movaps 6912(pB), rB2
            #elif KB == 146
               prefA((pfA))
            #endif

         #endif
         #if KB > 146
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6912(pA), rB0
            addps rB0, rC20
            #if KB > 147
               movaps 6928(pB), rB0
            #elif KB == 147
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6912(pA), rB1
            addps rB1, rC21
            #if KB > 147
               movaps 6944(pB), rB1
            #elif KB == 147
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 147
               movaps 6928(pA), rA0
            #elif KB == 147
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 147
               movaps 6944(pA), rA1
            #elif KB == 147
               prefA(-128(pfA))
            #endif

            mulps 6912(pA), rB2
            addps rB2, rC22
            #if KB > 147
               movaps 6960(pB), rB2
            #elif KB == 147
               prefA((pfA))
            #endif

         #endif
         #if KB > 147
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 6960(pA), rB0
            addps rB0, rC20
            #if KB > 148
               movaps 6976(pB), rB0
            #elif KB == 148
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 6960(pA), rB1
            addps rB1, rC21
            #if KB > 148
               movaps 6992(pB), rB1
            #elif KB == 148
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 148
               movaps 6976(pA), rA0
            #elif KB == 148
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 148
               movaps 6992(pA), rA1
            #elif KB == 148
               prefA(-128(pfA))
            #endif

            mulps 6960(pA), rB2
            addps rB2, rC22
            #if KB > 148
               movaps 7008(pB), rB2
            #elif KB == 148
               prefA((pfA))
            #endif

         #endif
         #if KB > 148
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7008(pA), rB0
            addps rB0, rC20
            #if KB > 149
               movaps 7024(pB), rB0
            #elif KB == 149
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7008(pA), rB1
            addps rB1, rC21
            #if KB > 149
               movaps 7040(pB), rB1
            #elif KB == 149
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 149
               movaps 7024(pA), rA0
            #elif KB == 149
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 149
               movaps 7040(pA), rA1
            #elif KB == 149
               prefA(-128(pfA))
            #endif

            mulps 7008(pA), rB2
            addps rB2, rC22
            #if KB > 149
               movaps 7056(pB), rB2
            #elif KB == 149
               prefA((pfA))
            #endif

         #endif
         #if KB > 149
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7056(pA), rB0
            addps rB0, rC20
            #if KB > 150
               movaps 7072(pB), rB0
            #elif KB == 150
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7056(pA), rB1
            addps rB1, rC21
            #if KB > 150
               movaps 7088(pB), rB1
            #elif KB == 150
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 150
               movaps 7072(pA), rA0
            #elif KB == 150
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 150
               movaps 7088(pA), rA1
            #elif KB == 150
               prefA(-128(pfA))
            #endif

            mulps 7056(pA), rB2
            addps rB2, rC22
            #if KB > 150
               movaps 7104(pB), rB2
            #elif KB == 150
               prefA((pfA))
            #endif

         #endif
         #if KB > 150
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7104(pA), rB0
            addps rB0, rC20
            #if KB > 151
               movaps 7120(pB), rB0
            #elif KB == 151
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7104(pA), rB1
            addps rB1, rC21
            #if KB > 151
               movaps 7136(pB), rB1
            #elif KB == 151
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 151
               movaps 7120(pA), rA0
            #elif KB == 151
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 151
               movaps 7136(pA), rA1
            #elif KB == 151
               prefA(-128(pfA))
            #endif

            mulps 7104(pA), rB2
            addps rB2, rC22
            #if KB > 151
               movaps 7152(pB), rB2
            #elif KB == 151
               prefA((pfA))
            #endif

         #endif
         #if KB > 151
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7152(pA), rB0
            addps rB0, rC20
            #if KB > 152
               movaps 7168(pB), rB0
            #elif KB == 152
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7152(pA), rB1
            addps rB1, rC21
            #if KB > 152
               movaps 7184(pB), rB1
            #elif KB == 152
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 152
               movaps 7168(pA), rA0
            #elif KB == 152
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 152
               movaps 7184(pA), rA1
            #elif KB == 152
               prefA(-128(pfA))
            #endif

            mulps 7152(pA), rB2
            addps rB2, rC22
            #if KB > 152
               movaps 7200(pB), rB2
            #elif KB == 152
               prefA((pfA))
            #endif

         #endif
         #if KB > 152
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7200(pA), rB0
            addps rB0, rC20
            #if KB > 153
               movaps 7216(pB), rB0
            #elif KB == 153
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7200(pA), rB1
            addps rB1, rC21
            #if KB > 153
               movaps 7232(pB), rB1
            #elif KB == 153
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 153
               movaps 7216(pA), rA0
            #elif KB == 153
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 153
               movaps 7232(pA), rA1
            #elif KB == 153
               prefA(-128(pfA))
            #endif

            mulps 7200(pA), rB2
            addps rB2, rC22
            #if KB > 153
               movaps 7248(pB), rB2
            #elif KB == 153
               prefA((pfA))
            #endif

         #endif
         #if KB > 153
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7248(pA), rB0
            addps rB0, rC20
            #if KB > 154
               movaps 7264(pB), rB0
            #elif KB == 154
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7248(pA), rB1
            addps rB1, rC21
            #if KB > 154
               movaps 7280(pB), rB1
            #elif KB == 154
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 154
               movaps 7264(pA), rA0
            #elif KB == 154
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 154
               movaps 7280(pA), rA1
            #elif KB == 154
               prefA(-128(pfA))
            #endif

            mulps 7248(pA), rB2
            addps rB2, rC22
            #if KB > 154
               movaps 7296(pB), rB2
            #elif KB == 154
               prefA((pfA))
            #endif

         #endif
         #if KB > 154
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7296(pA), rB0
            addps rB0, rC20
            #if KB > 155
               movaps 7312(pB), rB0
            #elif KB == 155
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7296(pA), rB1
            addps rB1, rC21
            #if KB > 155
               movaps 7328(pB), rB1
            #elif KB == 155
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 155
               movaps 7312(pA), rA0
            #elif KB == 155
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 155
               movaps 7328(pA), rA1
            #elif KB == 155
               prefA(-128(pfA))
            #endif

            mulps 7296(pA), rB2
            addps rB2, rC22
            #if KB > 155
               movaps 7344(pB), rB2
            #elif KB == 155
               prefA((pfA))
            #endif

         #endif
         #if KB > 155
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7344(pA), rB0
            addps rB0, rC20
            #if KB > 156
               movaps 7360(pB), rB0
            #elif KB == 156
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7344(pA), rB1
            addps rB1, rC21
            #if KB > 156
               movaps 7376(pB), rB1
            #elif KB == 156
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 156
               movaps 7360(pA), rA0
            #elif KB == 156
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 156
               movaps 7376(pA), rA1
            #elif KB == 156
               prefA(-128(pfA))
            #endif

            mulps 7344(pA), rB2
            addps rB2, rC22
            #if KB > 156
               movaps 7392(pB), rB2
            #elif KB == 156
               prefA((pfA))
            #endif

         #endif
         #if KB > 156
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7392(pA), rB0
            addps rB0, rC20
            #if KB > 157
               movaps 7408(pB), rB0
            #elif KB == 157
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7392(pA), rB1
            addps rB1, rC21
            #if KB > 157
               movaps 7424(pB), rB1
            #elif KB == 157
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 157
               movaps 7408(pA), rA0
            #elif KB == 157
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 157
               movaps 7424(pA), rA1
            #elif KB == 157
               prefA(-128(pfA))
            #endif

            mulps 7392(pA), rB2
            addps rB2, rC22
            #if KB > 157
               movaps 7440(pB), rB2
            #elif KB == 157
               prefA((pfA))
            #endif

         #endif
         #if KB > 157
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7440(pA), rB0
            addps rB0, rC20
            #if KB > 158
               movaps 7456(pB), rB0
            #elif KB == 158
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7440(pA), rB1
            addps rB1, rC21
            #if KB > 158
               movaps 7472(pB), rB1
            #elif KB == 158
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 158
               movaps 7456(pA), rA0
            #elif KB == 158
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 158
               movaps 7472(pA), rA1
            #elif KB == 158
               prefA(-128(pfA))
            #endif

            mulps 7440(pA), rB2
            addps rB2, rC22
            #if KB > 158
               movaps 7488(pB), rB2
            #elif KB == 158
               prefA((pfA))
            #endif

         #endif
         #if KB > 158
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7488(pA), rB0
            addps rB0, rC20
            #if KB > 159
               movaps 7504(pB), rB0
            #elif KB == 159
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7488(pA), rB1
            addps rB1, rC21
            #if KB > 159
               movaps 7520(pB), rB1
            #elif KB == 159
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 159
               movaps 7504(pA), rA0
            #elif KB == 159
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 159
               movaps 7520(pA), rA1
            #elif KB == 159
               prefA(-128(pfA))
            #endif

            mulps 7488(pA), rB2
            addps rB2, rC22
            #if KB > 159
               movaps 7536(pB), rB2
            #elif KB == 159
               prefA((pfA))
            #endif

         #endif
         #if KB > 159
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7536(pA), rB0
            addps rB0, rC20
            #if KB > 160
               movaps 7552(pB), rB0
            #elif KB == 160
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7536(pA), rB1
            addps rB1, rC21
            #if KB > 160
               movaps 7568(pB), rB1
            #elif KB == 160
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 160
               movaps 7552(pA), rA0
            #elif KB == 160
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 160
               movaps 7568(pA), rA1
            #elif KB == 160
               prefA(-128(pfA))
            #endif

            mulps 7536(pA), rB2
            addps rB2, rC22
            #if KB > 160
               movaps 7584(pB), rB2
            #elif KB == 160
               prefA((pfA))
            #endif

         #endif
         #if KB > 160
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7584(pA), rB0
            addps rB0, rC20
            #if KB > 161
               movaps 7600(pB), rB0
            #elif KB == 161
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7584(pA), rB1
            addps rB1, rC21
            #if KB > 161
               movaps 7616(pB), rB1
            #elif KB == 161
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 161
               movaps 7600(pA), rA0
            #elif KB == 161
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 161
               movaps 7616(pA), rA1
            #elif KB == 161
               prefA(-128(pfA))
            #endif

            mulps 7584(pA), rB2
            addps rB2, rC22
            #if KB > 161
               movaps 7632(pB), rB2
            #elif KB == 161
               prefA((pfA))
            #endif

         #endif
         #if KB > 161
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7632(pA), rB0
            addps rB0, rC20
            #if KB > 162
               movaps 7648(pB), rB0
            #elif KB == 162
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7632(pA), rB1
            addps rB1, rC21
            #if KB > 162
               movaps 7664(pB), rB1
            #elif KB == 162
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 162
               movaps 7648(pA), rA0
            #elif KB == 162
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 162
               movaps 7664(pA), rA1
            #elif KB == 162
               prefA(-128(pfA))
            #endif

            mulps 7632(pA), rB2
            addps rB2, rC22
            #if KB > 162
               movaps 7680(pB), rB2
            #elif KB == 162
               prefA((pfA))
            #endif

         #endif
         #if KB > 162
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7680(pA), rB0
            addps rB0, rC20
            #if KB > 163
               movaps 7696(pB), rB0
            #elif KB == 163
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7680(pA), rB1
            addps rB1, rC21
            #if KB > 163
               movaps 7712(pB), rB1
            #elif KB == 163
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 163
               movaps 7696(pA), rA0
            #elif KB == 163
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 163
               movaps 7712(pA), rA1
            #elif KB == 163
               prefA(-128(pfA))
            #endif

            mulps 7680(pA), rB2
            addps rB2, rC22
            #if KB > 163
               movaps 7728(pB), rB2
            #elif KB == 163
               prefA((pfA))
            #endif

         #endif
         #if KB > 163
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7728(pA), rB0
            addps rB0, rC20
            #if KB > 164
               movaps 7744(pB), rB0
            #elif KB == 164
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7728(pA), rB1
            addps rB1, rC21
            #if KB > 164
               movaps 7760(pB), rB1
            #elif KB == 164
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 164
               movaps 7744(pA), rA0
            #elif KB == 164
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 164
               movaps 7760(pA), rA1
            #elif KB == 164
               prefA(-128(pfA))
            #endif

            mulps 7728(pA), rB2
            addps rB2, rC22
            #if KB > 164
               movaps 7776(pB), rB2
            #elif KB == 164
               prefA((pfA))
            #endif

         #endif
         #if KB > 164
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7776(pA), rB0
            addps rB0, rC20
            #if KB > 165
               movaps 7792(pB), rB0
            #elif KB == 165
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7776(pA), rB1
            addps rB1, rC21
            #if KB > 165
               movaps 7808(pB), rB1
            #elif KB == 165
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 165
               movaps 7792(pA), rA0
            #elif KB == 165
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 165
               movaps 7808(pA), rA1
            #elif KB == 165
               prefA(-128(pfA))
            #endif

            mulps 7776(pA), rB2
            addps rB2, rC22
            #if KB > 165
               movaps 7824(pB), rB2
            #elif KB == 165
               prefA((pfA))
            #endif

         #endif
         #if KB > 165
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7824(pA), rB0
            addps rB0, rC20
            #if KB > 166
               movaps 7840(pB), rB0
            #elif KB == 166
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7824(pA), rB1
            addps rB1, rC21
            #if KB > 166
               movaps 7856(pB), rB1
            #elif KB == 166
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 166
               movaps 7840(pA), rA0
            #elif KB == 166
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 166
               movaps 7856(pA), rA1
            #elif KB == 166
               prefA(-128(pfA))
            #endif

            mulps 7824(pA), rB2
            addps rB2, rC22
            #if KB > 166
               movaps 7872(pB), rB2
            #elif KB == 166
               prefA((pfA))
            #endif

         #endif
         #if KB > 166
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7872(pA), rB0
            addps rB0, rC20
            #if KB > 167
               movaps 7888(pB), rB0
            #elif KB == 167
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7872(pA), rB1
            addps rB1, rC21
            #if KB > 167
               movaps 7904(pB), rB1
            #elif KB == 167
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 167
               movaps 7888(pA), rA0
            #elif KB == 167
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 167
               movaps 7904(pA), rA1
            #elif KB == 167
               prefA(-128(pfA))
            #endif

            mulps 7872(pA), rB2
            addps rB2, rC22
            #if KB > 167
               movaps 7920(pB), rB2
            #elif KB == 167
               prefA((pfA))
            #endif

         #endif
         #if KB > 167
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7920(pA), rB0
            addps rB0, rC20
            #if KB > 168
               movaps 7936(pB), rB0
            #elif KB == 168
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7920(pA), rB1
            addps rB1, rC21
            #if KB > 168
               movaps 7952(pB), rB1
            #elif KB == 168
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 168
               movaps 7936(pA), rA0
            #elif KB == 168
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 168
               movaps 7952(pA), rA1
            #elif KB == 168
               prefA(-128(pfA))
            #endif

            mulps 7920(pA), rB2
            addps rB2, rC22
            #if KB > 168
               movaps 7968(pB), rB2
            #elif KB == 168
               prefA((pfA))
            #endif

         #endif
         #if KB > 168
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 7968(pA), rB0
            addps rB0, rC20
            #if KB > 169
               movaps 7984(pB), rB0
            #elif KB == 169
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 7968(pA), rB1
            addps rB1, rC21
            #if KB > 169
               movaps 8000(pB), rB1
            #elif KB == 169
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 169
               movaps 7984(pA), rA0
            #elif KB == 169
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 169
               movaps 8000(pA), rA1
            #elif KB == 169
               prefA(-128(pfA))
            #endif

            mulps 7968(pA), rB2
            addps rB2, rC22
            #if KB > 169
               movaps 8016(pB), rB2
            #elif KB == 169
               prefA((pfA))
            #endif

         #endif
         #if KB > 169
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8016(pA), rB0
            addps rB0, rC20
            #if KB > 170
               movaps 8032(pB), rB0
            #elif KB == 170
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8016(pA), rB1
            addps rB1, rC21
            #if KB > 170
               movaps 8048(pB), rB1
            #elif KB == 170
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 170
               movaps 8032(pA), rA0
            #elif KB == 170
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 170
               movaps 8048(pA), rA1
            #elif KB == 170
               prefA(-128(pfA))
            #endif

            mulps 8016(pA), rB2
            addps rB2, rC22
            #if KB > 170
               movaps 8064(pB), rB2
            #elif KB == 170
               prefA((pfA))
            #endif

         #endif
         #if KB > 170
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8064(pA), rB0
            addps rB0, rC20
            #if KB > 171
               movaps 8080(pB), rB0
            #elif KB == 171
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8064(pA), rB1
            addps rB1, rC21
            #if KB > 171
               movaps 8096(pB), rB1
            #elif KB == 171
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 171
               movaps 8080(pA), rA0
            #elif KB == 171
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 171
               movaps 8096(pA), rA1
            #elif KB == 171
               prefA(-128(pfA))
            #endif

            mulps 8064(pA), rB2
            addps rB2, rC22
            #if KB > 171
               movaps 8112(pB), rB2
            #elif KB == 171
               prefA((pfA))
            #endif

         #endif
         #if KB > 171
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8112(pA), rB0
            addps rB0, rC20
            #if KB > 172
               movaps 8128(pB), rB0
            #elif KB == 172
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8112(pA), rB1
            addps rB1, rC21
            #if KB > 172
               movaps 8144(pB), rB1
            #elif KB == 172
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 172
               movaps 8128(pA), rA0
            #elif KB == 172
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 172
               movaps 8144(pA), rA1
            #elif KB == 172
               prefA(-128(pfA))
            #endif

            mulps 8112(pA), rB2
            addps rB2, rC22
            #if KB > 172
               movaps 8160(pB), rB2
            #elif KB == 172
               prefA((pfA))
            #endif

         #endif
         #if KB > 172
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8160(pA), rB0
            addps rB0, rC20
            #if KB > 173
               movaps 8176(pB), rB0
            #elif KB == 173
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8160(pA), rB1
            addps rB1, rC21
            #if KB > 173
               movaps 8192(pB), rB1
            #elif KB == 173
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 173
               movaps 8176(pA), rA0
            #elif KB == 173
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 173
               movaps 8192(pA), rA1
            #elif KB == 173
               prefA(-128(pfA))
            #endif

            mulps 8160(pA), rB2
            addps rB2, rC22
            #if KB > 173
               movaps 8208(pB), rB2
            #elif KB == 173
               prefA((pfA))
            #endif

         #endif
         #if KB > 173
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8208(pA), rB0
            addps rB0, rC20
            #if KB > 174
               movaps 8224(pB), rB0
            #elif KB == 174
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8208(pA), rB1
            addps rB1, rC21
            #if KB > 174
               movaps 8240(pB), rB1
            #elif KB == 174
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 174
               movaps 8224(pA), rA0
            #elif KB == 174
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 174
               movaps 8240(pA), rA1
            #elif KB == 174
               prefA(-128(pfA))
            #endif

            mulps 8208(pA), rB2
            addps rB2, rC22
            #if KB > 174
               movaps 8256(pB), rB2
            #elif KB == 174
               prefA((pfA))
            #endif

         #endif
         #if KB > 174
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8256(pA), rB0
            addps rB0, rC20
            #if KB > 175
               movaps 8272(pB), rB0
            #elif KB == 175
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8256(pA), rB1
            addps rB1, rC21
            #if KB > 175
               movaps 8288(pB), rB1
            #elif KB == 175
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 175
               movaps 8272(pA), rA0
            #elif KB == 175
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 175
               movaps 8288(pA), rA1
            #elif KB == 175
               prefA(-128(pfA))
            #endif

            mulps 8256(pA), rB2
            addps rB2, rC22
            #if KB > 175
               movaps 8304(pB), rB2
            #elif KB == 175
               prefA((pfA))
            #endif

         #endif
         #if KB > 175
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8304(pA), rB0
            addps rB0, rC20
            #if KB > 176
               movaps 8320(pB), rB0
            #elif KB == 176
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8304(pA), rB1
            addps rB1, rC21
            #if KB > 176
               movaps 8336(pB), rB1
            #elif KB == 176
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 176
               movaps 8320(pA), rA0
            #elif KB == 176
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 176
               movaps 8336(pA), rA1
            #elif KB == 176
               prefA(-128(pfA))
            #endif

            mulps 8304(pA), rB2
            addps rB2, rC22
            #if KB > 176
               movaps 8352(pB), rB2
            #elif KB == 176
               prefA((pfA))
            #endif

         #endif
         #if KB > 176
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8352(pA), rB0
            addps rB0, rC20
            #if KB > 177
               movaps 8368(pB), rB0
            #elif KB == 177
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8352(pA), rB1
            addps rB1, rC21
            #if KB > 177
               movaps 8384(pB), rB1
            #elif KB == 177
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 177
               movaps 8368(pA), rA0
            #elif KB == 177
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 177
               movaps 8384(pA), rA1
            #elif KB == 177
               prefA(-128(pfA))
            #endif

            mulps 8352(pA), rB2
            addps rB2, rC22
            #if KB > 177
               movaps 8400(pB), rB2
            #elif KB == 177
               prefA((pfA))
            #endif

         #endif
         #if KB > 177
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8400(pA), rB0
            addps rB0, rC20
            #if KB > 178
               movaps 8416(pB), rB0
            #elif KB == 178
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8400(pA), rB1
            addps rB1, rC21
            #if KB > 178
               movaps 8432(pB), rB1
            #elif KB == 178
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 178
               movaps 8416(pA), rA0
            #elif KB == 178
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 178
               movaps 8432(pA), rA1
            #elif KB == 178
               prefA(-128(pfA))
            #endif

            mulps 8400(pA), rB2
            addps rB2, rC22
            #if KB > 178
               movaps 8448(pB), rB2
            #elif KB == 178
               prefA((pfA))
            #endif

         #endif
         #if KB > 178
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8448(pA), rB0
            addps rB0, rC20
            #if KB > 179
               movaps 8464(pB), rB0
            #elif KB == 179
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8448(pA), rB1
            addps rB1, rC21
            #if KB > 179
               movaps 8480(pB), rB1
            #elif KB == 179
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 179
               movaps 8464(pA), rA0
            #elif KB == 179
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 179
               movaps 8480(pA), rA1
            #elif KB == 179
               prefA(-128(pfA))
            #endif

            mulps 8448(pA), rB2
            addps rB2, rC22
            #if KB > 179
               movaps 8496(pB), rB2
            #elif KB == 179
               prefA((pfA))
            #endif

         #endif
         #if KB > 179
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8496(pA), rB0
            addps rB0, rC20
            #if KB > 180
               movaps 8512(pB), rB0
            #elif KB == 180
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8496(pA), rB1
            addps rB1, rC21
            #if KB > 180
               movaps 8528(pB), rB1
            #elif KB == 180
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 180
               movaps 8512(pA), rA0
            #elif KB == 180
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 180
               movaps 8528(pA), rA1
            #elif KB == 180
               prefA(-128(pfA))
            #endif

            mulps 8496(pA), rB2
            addps rB2, rC22
            #if KB > 180
               movaps 8544(pB), rB2
            #elif KB == 180
               prefA((pfA))
            #endif

         #endif
         #if KB > 180
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8544(pA), rB0
            addps rB0, rC20
            #if KB > 181
               movaps 8560(pB), rB0
            #elif KB == 181
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8544(pA), rB1
            addps rB1, rC21
            #if KB > 181
               movaps 8576(pB), rB1
            #elif KB == 181
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 181
               movaps 8560(pA), rA0
            #elif KB == 181
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 181
               movaps 8576(pA), rA1
            #elif KB == 181
               prefA(-128(pfA))
            #endif

            mulps 8544(pA), rB2
            addps rB2, rC22
            #if KB > 181
               movaps 8592(pB), rB2
            #elif KB == 181
               prefA((pfA))
            #endif

         #endif
         #if KB > 181
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8592(pA), rB0
            addps rB0, rC20
            #if KB > 182
               movaps 8608(pB), rB0
            #elif KB == 182
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8592(pA), rB1
            addps rB1, rC21
            #if KB > 182
               movaps 8624(pB), rB1
            #elif KB == 182
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 182
               movaps 8608(pA), rA0
            #elif KB == 182
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 182
               movaps 8624(pA), rA1
            #elif KB == 182
               prefA(-128(pfA))
            #endif

            mulps 8592(pA), rB2
            addps rB2, rC22
            #if KB > 182
               movaps 8640(pB), rB2
            #elif KB == 182
               prefA((pfA))
            #endif

         #endif
         #if KB > 182
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8640(pA), rB0
            addps rB0, rC20
            #if KB > 183
               movaps 8656(pB), rB0
            #elif KB == 183
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8640(pA), rB1
            addps rB1, rC21
            #if KB > 183
               movaps 8672(pB), rB1
            #elif KB == 183
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 183
               movaps 8656(pA), rA0
            #elif KB == 183
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 183
               movaps 8672(pA), rA1
            #elif KB == 183
               prefA(-128(pfA))
            #endif

            mulps 8640(pA), rB2
            addps rB2, rC22
            #if KB > 183
               movaps 8688(pB), rB2
            #elif KB == 183
               prefA((pfA))
            #endif

         #endif
         #if KB > 183
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8688(pA), rB0
            addps rB0, rC20
            #if KB > 184
               movaps 8704(pB), rB0
            #elif KB == 184
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8688(pA), rB1
            addps rB1, rC21
            #if KB > 184
               movaps 8720(pB), rB1
            #elif KB == 184
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 184
               movaps 8704(pA), rA0
            #elif KB == 184
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 184
               movaps 8720(pA), rA1
            #elif KB == 184
               prefA(-128(pfA))
            #endif

            mulps 8688(pA), rB2
            addps rB2, rC22
            #if KB > 184
               movaps 8736(pB), rB2
            #elif KB == 184
               prefA((pfA))
            #endif

         #endif
         #if KB > 184
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8736(pA), rB0
            addps rB0, rC20
            #if KB > 185
               movaps 8752(pB), rB0
            #elif KB == 185
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8736(pA), rB1
            addps rB1, rC21
            #if KB > 185
               movaps 8768(pB), rB1
            #elif KB == 185
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 185
               movaps 8752(pA), rA0
            #elif KB == 185
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 185
               movaps 8768(pA), rA1
            #elif KB == 185
               prefA(-128(pfA))
            #endif

            mulps 8736(pA), rB2
            addps rB2, rC22
            #if KB > 185
               movaps 8784(pB), rB2
            #elif KB == 185
               prefA((pfA))
            #endif

         #endif
         #if KB > 185
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8784(pA), rB0
            addps rB0, rC20
            #if KB > 186
               movaps 8800(pB), rB0
            #elif KB == 186
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8784(pA), rB1
            addps rB1, rC21
            #if KB > 186
               movaps 8816(pB), rB1
            #elif KB == 186
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 186
               movaps 8800(pA), rA0
            #elif KB == 186
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 186
               movaps 8816(pA), rA1
            #elif KB == 186
               prefA(-128(pfA))
            #endif

            mulps 8784(pA), rB2
            addps rB2, rC22
            #if KB > 186
               movaps 8832(pB), rB2
            #elif KB == 186
               prefA((pfA))
            #endif

         #endif
         #if KB > 186
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8832(pA), rB0
            addps rB0, rC20
            #if KB > 187
               movaps 8848(pB), rB0
            #elif KB == 187
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8832(pA), rB1
            addps rB1, rC21
            #if KB > 187
               movaps 8864(pB), rB1
            #elif KB == 187
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 187
               movaps 8848(pA), rA0
            #elif KB == 187
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 187
               movaps 8864(pA), rA1
            #elif KB == 187
               prefA(-128(pfA))
            #endif

            mulps 8832(pA), rB2
            addps rB2, rC22
            #if KB > 187
               movaps 8880(pB), rB2
            #elif KB == 187
               prefA((pfA))
            #endif

         #endif
         #if KB > 187
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8880(pA), rB0
            addps rB0, rC20
            #if KB > 188
               movaps 8896(pB), rB0
            #elif KB == 188
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8880(pA), rB1
            addps rB1, rC21
            #if KB > 188
               movaps 8912(pB), rB1
            #elif KB == 188
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 188
               movaps 8896(pA), rA0
            #elif KB == 188
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 188
               movaps 8912(pA), rA1
            #elif KB == 188
               prefA(-128(pfA))
            #endif

            mulps 8880(pA), rB2
            addps rB2, rC22
            #if KB > 188
               movaps 8928(pB), rB2
            #elif KB == 188
               prefA((pfA))
            #endif

         #endif
         #if KB > 188
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8928(pA), rB0
            addps rB0, rC20
            #if KB > 189
               movaps 8944(pB), rB0
            #elif KB == 189
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8928(pA), rB1
            addps rB1, rC21
            #if KB > 189
               movaps 8960(pB), rB1
            #elif KB == 189
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 189
               movaps 8944(pA), rA0
            #elif KB == 189
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 189
               movaps 8960(pA), rA1
            #elif KB == 189
               prefA(-128(pfA))
            #endif

            mulps 8928(pA), rB2
            addps rB2, rC22
            #if KB > 189
               movaps 8976(pB), rB2
            #elif KB == 189
               prefA((pfA))
            #endif

         #endif
         #if KB > 189
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 8976(pA), rB0
            addps rB0, rC20
            #if KB > 190
               movaps 8992(pB), rB0
            #elif KB == 190
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 8976(pA), rB1
            addps rB1, rC21
            #if KB > 190
               movaps 9008(pB), rB1
            #elif KB == 190
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 190
               movaps 8992(pA), rA0
            #elif KB == 190
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 190
               movaps 9008(pA), rA1
            #elif KB == 190
               prefA(-128(pfA))
            #endif

            mulps 8976(pA), rB2
            addps rB2, rC22
            #if KB > 190
               movaps 9024(pB), rB2
            #elif KB == 190
               prefA((pfA))
            #endif

         #endif
         #if KB > 190
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9024(pA), rB0
            addps rB0, rC20
            #if KB > 191
               movaps 9040(pB), rB0
            #elif KB == 191
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9024(pA), rB1
            addps rB1, rC21
            #if KB > 191
               movaps 9056(pB), rB1
            #elif KB == 191
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 191
               movaps 9040(pA), rA0
            #elif KB == 191
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 191
               movaps 9056(pA), rA1
            #elif KB == 191
               prefA(-128(pfA))
            #endif

            mulps 9024(pA), rB2
            addps rB2, rC22
            #if KB > 191
               movaps 9072(pB), rB2
            #elif KB == 191
               prefA((pfA))
            #endif

         #endif
         #if KB > 191
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9072(pA), rB0
            addps rB0, rC20
            #if KB > 192
               movaps 9088(pB), rB0
            #elif KB == 192
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9072(pA), rB1
            addps rB1, rC21
            #if KB > 192
               movaps 9104(pB), rB1
            #elif KB == 192
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 192
               movaps 9088(pA), rA0
            #elif KB == 192
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 192
               movaps 9104(pA), rA1
            #elif KB == 192
               prefA(-128(pfA))
            #endif

            mulps 9072(pA), rB2
            addps rB2, rC22
            #if KB > 192
               movaps 9120(pB), rB2
            #elif KB == 192
               prefA((pfA))
            #endif

         #endif
         #if KB > 192
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9120(pA), rB0
            addps rB0, rC20
            #if KB > 193
               movaps 9136(pB), rB0
            #elif KB == 193
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9120(pA), rB1
            addps rB1, rC21
            #if KB > 193
               movaps 9152(pB), rB1
            #elif KB == 193
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 193
               movaps 9136(pA), rA0
            #elif KB == 193
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 193
               movaps 9152(pA), rA1
            #elif KB == 193
               prefA(-128(pfA))
            #endif

            mulps 9120(pA), rB2
            addps rB2, rC22
            #if KB > 193
               movaps 9168(pB), rB2
            #elif KB == 193
               prefA((pfA))
            #endif

         #endif
         #if KB > 193
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9168(pA), rB0
            addps rB0, rC20
            #if KB > 194
               movaps 9184(pB), rB0
            #elif KB == 194
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9168(pA), rB1
            addps rB1, rC21
            #if KB > 194
               movaps 9200(pB), rB1
            #elif KB == 194
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 194
               movaps 9184(pA), rA0
            #elif KB == 194
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 194
               movaps 9200(pA), rA1
            #elif KB == 194
               prefA(-128(pfA))
            #endif

            mulps 9168(pA), rB2
            addps rB2, rC22
            #if KB > 194
               movaps 9216(pB), rB2
            #elif KB == 194
               prefA((pfA))
            #endif

         #endif
         #if KB > 194
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9216(pA), rB0
            addps rB0, rC20
            #if KB > 195
               movaps 9232(pB), rB0
            #elif KB == 195
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9216(pA), rB1
            addps rB1, rC21
            #if KB > 195
               movaps 9248(pB), rB1
            #elif KB == 195
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 195
               movaps 9232(pA), rA0
            #elif KB == 195
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 195
               movaps 9248(pA), rA1
            #elif KB == 195
               prefA(-128(pfA))
            #endif

            mulps 9216(pA), rB2
            addps rB2, rC22
            #if KB > 195
               movaps 9264(pB), rB2
            #elif KB == 195
               prefA((pfA))
            #endif

         #endif
         #if KB > 195
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9264(pA), rB0
            addps rB0, rC20
            #if KB > 196
               movaps 9280(pB), rB0
            #elif KB == 196
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9264(pA), rB1
            addps rB1, rC21
            #if KB > 196
               movaps 9296(pB), rB1
            #elif KB == 196
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 196
               movaps 9280(pA), rA0
            #elif KB == 196
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 196
               movaps 9296(pA), rA1
            #elif KB == 196
               prefA(-128(pfA))
            #endif

            mulps 9264(pA), rB2
            addps rB2, rC22
            #if KB > 196
               movaps 9312(pB), rB2
            #elif KB == 196
               prefA((pfA))
            #endif

         #endif
         #if KB > 196
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9312(pA), rB0
            addps rB0, rC20
            #if KB > 197
               movaps 9328(pB), rB0
            #elif KB == 197
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9312(pA), rB1
            addps rB1, rC21
            #if KB > 197
               movaps 9344(pB), rB1
            #elif KB == 197
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 197
               movaps 9328(pA), rA0
            #elif KB == 197
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 197
               movaps 9344(pA), rA1
            #elif KB == 197
               prefA(-128(pfA))
            #endif

            mulps 9312(pA), rB2
            addps rB2, rC22
            #if KB > 197
               movaps 9360(pB), rB2
            #elif KB == 197
               prefA((pfA))
            #endif

         #endif
         #if KB > 197
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9360(pA), rB0
            addps rB0, rC20
            #if KB > 198
               movaps 9376(pB), rB0
            #elif KB == 198
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9360(pA), rB1
            addps rB1, rC21
            #if KB > 198
               movaps 9392(pB), rB1
            #elif KB == 198
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 198
               movaps 9376(pA), rA0
            #elif KB == 198
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 198
               movaps 9392(pA), rA1
            #elif KB == 198
               prefA(-128(pfA))
            #endif

            mulps 9360(pA), rB2
            addps rB2, rC22
            #if KB > 198
               movaps 9408(pB), rB2
            #elif KB == 198
               prefA((pfA))
            #endif

         #endif
         #if KB > 198
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9408(pA), rB0
            addps rB0, rC20
            #if KB > 199
               movaps 9424(pB), rB0
            #elif KB == 199
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9408(pA), rB1
            addps rB1, rC21
            #if KB > 199
               movaps 9440(pB), rB1
            #elif KB == 199
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 199
               movaps 9424(pA), rA0
            #elif KB == 199
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 199
               movaps 9440(pA), rA1
            #elif KB == 199
               prefA(-128(pfA))
            #endif

            mulps 9408(pA), rB2
            addps rB2, rC22
            #if KB > 199
               movaps 9456(pB), rB2
            #elif KB == 199
               prefA((pfA))
            #endif

         #endif
         #if KB > 199
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9456(pA), rB0
            addps rB0, rC20
            #if KB > 200
               movaps 9472(pB), rB0
            #elif KB == 200
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9456(pA), rB1
            addps rB1, rC21
            #if KB > 200
               movaps 9488(pB), rB1
            #elif KB == 200
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 200
               movaps 9472(pA), rA0
            #elif KB == 200
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 200
               movaps 9488(pA), rA1
            #elif KB == 200
               prefA(-128(pfA))
            #endif

            mulps 9456(pA), rB2
            addps rB2, rC22
            #if KB > 200
               movaps 9504(pB), rB2
            #elif KB == 200
               prefA((pfA))
            #endif

         #endif
         #if KB > 200
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9504(pA), rB0
            addps rB0, rC20
            #if KB > 201
               movaps 9520(pB), rB0
            #elif KB == 201
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9504(pA), rB1
            addps rB1, rC21
            #if KB > 201
               movaps 9536(pB), rB1
            #elif KB == 201
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 201
               movaps 9520(pA), rA0
            #elif KB == 201
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 201
               movaps 9536(pA), rA1
            #elif KB == 201
               prefA(-128(pfA))
            #endif

            mulps 9504(pA), rB2
            addps rB2, rC22
            #if KB > 201
               movaps 9552(pB), rB2
            #elif KB == 201
               prefA((pfA))
            #endif

         #endif
         #if KB > 201
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9552(pA), rB0
            addps rB0, rC20
            #if KB > 202
               movaps 9568(pB), rB0
            #elif KB == 202
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9552(pA), rB1
            addps rB1, rC21
            #if KB > 202
               movaps 9584(pB), rB1
            #elif KB == 202
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 202
               movaps 9568(pA), rA0
            #elif KB == 202
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 202
               movaps 9584(pA), rA1
            #elif KB == 202
               prefA(-128(pfA))
            #endif

            mulps 9552(pA), rB2
            addps rB2, rC22
            #if KB > 202
               movaps 9600(pB), rB2
            #elif KB == 202
               prefA((pfA))
            #endif

         #endif
         #if KB > 202
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9600(pA), rB0
            addps rB0, rC20
            #if KB > 203
               movaps 9616(pB), rB0
            #elif KB == 203
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9600(pA), rB1
            addps rB1, rC21
            #if KB > 203
               movaps 9632(pB), rB1
            #elif KB == 203
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 203
               movaps 9616(pA), rA0
            #elif KB == 203
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 203
               movaps 9632(pA), rA1
            #elif KB == 203
               prefA(-128(pfA))
            #endif

            mulps 9600(pA), rB2
            addps rB2, rC22
            #if KB > 203
               movaps 9648(pB), rB2
            #elif KB == 203
               prefA((pfA))
            #endif

         #endif
         #if KB > 203
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9648(pA), rB0
            addps rB0, rC20
            #if KB > 204
               movaps 9664(pB), rB0
            #elif KB == 204
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9648(pA), rB1
            addps rB1, rC21
            #if KB > 204
               movaps 9680(pB), rB1
            #elif KB == 204
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 204
               movaps 9664(pA), rA0
            #elif KB == 204
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 204
               movaps 9680(pA), rA1
            #elif KB == 204
               prefA(-128(pfA))
            #endif

            mulps 9648(pA), rB2
            addps rB2, rC22
            #if KB > 204
               movaps 9696(pB), rB2
            #elif KB == 204
               prefA((pfA))
            #endif

         #endif
         #if KB > 204
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9696(pA), rB0
            addps rB0, rC20
            #if KB > 205
               movaps 9712(pB), rB0
            #elif KB == 205
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9696(pA), rB1
            addps rB1, rC21
            #if KB > 205
               movaps 9728(pB), rB1
            #elif KB == 205
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 205
               movaps 9712(pA), rA0
            #elif KB == 205
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 205
               movaps 9728(pA), rA1
            #elif KB == 205
               prefA(-128(pfA))
            #endif

            mulps 9696(pA), rB2
            addps rB2, rC22
            #if KB > 205
               movaps 9744(pB), rB2
            #elif KB == 205
               prefA((pfA))
            #endif

         #endif
         #if KB > 205
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9744(pA), rB0
            addps rB0, rC20
            #if KB > 206
               movaps 9760(pB), rB0
            #elif KB == 206
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9744(pA), rB1
            addps rB1, rC21
            #if KB > 206
               movaps 9776(pB), rB1
            #elif KB == 206
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 206
               movaps 9760(pA), rA0
            #elif KB == 206
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 206
               movaps 9776(pA), rA1
            #elif KB == 206
               prefA(-128(pfA))
            #endif

            mulps 9744(pA), rB2
            addps rB2, rC22
            #if KB > 206
               movaps 9792(pB), rB2
            #elif KB == 206
               prefA((pfA))
            #endif

         #endif
         #if KB > 206
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9792(pA), rB0
            addps rB0, rC20
            #if KB > 207
               movaps 9808(pB), rB0
            #elif KB == 207
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9792(pA), rB1
            addps rB1, rC21
            #if KB > 207
               movaps 9824(pB), rB1
            #elif KB == 207
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 207
               movaps 9808(pA), rA0
            #elif KB == 207
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 207
               movaps 9824(pA), rA1
            #elif KB == 207
               prefA(-128(pfA))
            #endif

            mulps 9792(pA), rB2
            addps rB2, rC22
            #if KB > 207
               movaps 9840(pB), rB2
            #elif KB == 207
               prefA((pfA))
            #endif

         #endif
         #if KB > 207
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9840(pA), rB0
            addps rB0, rC20
            #if KB > 208
               movaps 9856(pB), rB0
            #elif KB == 208
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9840(pA), rB1
            addps rB1, rC21
            #if KB > 208
               movaps 9872(pB), rB1
            #elif KB == 208
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 208
               movaps 9856(pA), rA0
            #elif KB == 208
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 208
               movaps 9872(pA), rA1
            #elif KB == 208
               prefA(-128(pfA))
            #endif

            mulps 9840(pA), rB2
            addps rB2, rC22
            #if KB > 208
               movaps 9888(pB), rB2
            #elif KB == 208
               prefA((pfA))
            #endif

         #endif
         #if KB > 208
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9888(pA), rB0
            addps rB0, rC20
            #if KB > 209
               movaps 9904(pB), rB0
            #elif KB == 209
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9888(pA), rB1
            addps rB1, rC21
            #if KB > 209
               movaps 9920(pB), rB1
            #elif KB == 209
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 209
               movaps 9904(pA), rA0
            #elif KB == 209
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 209
               movaps 9920(pA), rA1
            #elif KB == 209
               prefA(-128(pfA))
            #endif

            mulps 9888(pA), rB2
            addps rB2, rC22
            #if KB > 209
               movaps 9936(pB), rB2
            #elif KB == 209
               prefA((pfA))
            #endif

         #endif
         #if KB > 209
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9936(pA), rB0
            addps rB0, rC20
            #if KB > 210
               movaps 9952(pB), rB0
            #elif KB == 210
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9936(pA), rB1
            addps rB1, rC21
            #if KB > 210
               movaps 9968(pB), rB1
            #elif KB == 210
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 210
               movaps 9952(pA), rA0
            #elif KB == 210
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 210
               movaps 9968(pA), rA1
            #elif KB == 210
               prefA(-128(pfA))
            #endif

            mulps 9936(pA), rB2
            addps rB2, rC22
            #if KB > 210
               movaps 9984(pB), rB2
            #elif KB == 210
               prefA((pfA))
            #endif

         #endif
         #if KB > 210
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 9984(pA), rB0
            addps rB0, rC20
            #if KB > 211
               movaps 10000(pB), rB0
            #elif KB == 211
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 9984(pA), rB1
            addps rB1, rC21
            #if KB > 211
               movaps 10016(pB), rB1
            #elif KB == 211
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 211
               movaps 10000(pA), rA0
            #elif KB == 211
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 211
               movaps 10016(pA), rA1
            #elif KB == 211
               prefA(-128(pfA))
            #endif

            mulps 9984(pA), rB2
            addps rB2, rC22
            #if KB > 211
               movaps 10032(pB), rB2
            #elif KB == 211
               prefA((pfA))
            #endif

         #endif
         #if KB > 211
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10032(pA), rB0
            addps rB0, rC20
            #if KB > 212
               movaps 10048(pB), rB0
            #elif KB == 212
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10032(pA), rB1
            addps rB1, rC21
            #if KB > 212
               movaps 10064(pB), rB1
            #elif KB == 212
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 212
               movaps 10048(pA), rA0
            #elif KB == 212
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 212
               movaps 10064(pA), rA1
            #elif KB == 212
               prefA(-128(pfA))
            #endif

            mulps 10032(pA), rB2
            addps rB2, rC22
            #if KB > 212
               movaps 10080(pB), rB2
            #elif KB == 212
               prefA((pfA))
            #endif

         #endif
         #if KB > 212
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10080(pA), rB0
            addps rB0, rC20
            #if KB > 213
               movaps 10096(pB), rB0
            #elif KB == 213
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10080(pA), rB1
            addps rB1, rC21
            #if KB > 213
               movaps 10112(pB), rB1
            #elif KB == 213
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 213
               movaps 10096(pA), rA0
            #elif KB == 213
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 213
               movaps 10112(pA), rA1
            #elif KB == 213
               prefA(-128(pfA))
            #endif

            mulps 10080(pA), rB2
            addps rB2, rC22
            #if KB > 213
               movaps 10128(pB), rB2
            #elif KB == 213
               prefA((pfA))
            #endif

         #endif
         #if KB > 213
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10128(pA), rB0
            addps rB0, rC20
            #if KB > 214
               movaps 10144(pB), rB0
            #elif KB == 214
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10128(pA), rB1
            addps rB1, rC21
            #if KB > 214
               movaps 10160(pB), rB1
            #elif KB == 214
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 214
               movaps 10144(pA), rA0
            #elif KB == 214
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 214
               movaps 10160(pA), rA1
            #elif KB == 214
               prefA(-128(pfA))
            #endif

            mulps 10128(pA), rB2
            addps rB2, rC22
            #if KB > 214
               movaps 10176(pB), rB2
            #elif KB == 214
               prefA((pfA))
            #endif

         #endif
         #if KB > 214
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10176(pA), rB0
            addps rB0, rC20
            #if KB > 215
               movaps 10192(pB), rB0
            #elif KB == 215
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10176(pA), rB1
            addps rB1, rC21
            #if KB > 215
               movaps 10208(pB), rB1
            #elif KB == 215
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 215
               movaps 10192(pA), rA0
            #elif KB == 215
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 215
               movaps 10208(pA), rA1
            #elif KB == 215
               prefA(-128(pfA))
            #endif

            mulps 10176(pA), rB2
            addps rB2, rC22
            #if KB > 215
               movaps 10224(pB), rB2
            #elif KB == 215
               prefA((pfA))
            #endif

         #endif
         #if KB > 215
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10224(pA), rB0
            addps rB0, rC20
            #if KB > 216
               movaps 10240(pB), rB0
            #elif KB == 216
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10224(pA), rB1
            addps rB1, rC21
            #if KB > 216
               movaps 10256(pB), rB1
            #elif KB == 216
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 216
               movaps 10240(pA), rA0
            #elif KB == 216
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 216
               movaps 10256(pA), rA1
            #elif KB == 216
               prefA(-128(pfA))
            #endif

            mulps 10224(pA), rB2
            addps rB2, rC22
            #if KB > 216
               movaps 10272(pB), rB2
            #elif KB == 216
               prefA((pfA))
            #endif

         #endif
         #if KB > 216
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10272(pA), rB0
            addps rB0, rC20
            #if KB > 217
               movaps 10288(pB), rB0
            #elif KB == 217
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10272(pA), rB1
            addps rB1, rC21
            #if KB > 217
               movaps 10304(pB), rB1
            #elif KB == 217
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 217
               movaps 10288(pA), rA0
            #elif KB == 217
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 217
               movaps 10304(pA), rA1
            #elif KB == 217
               prefA(-128(pfA))
            #endif

            mulps 10272(pA), rB2
            addps rB2, rC22
            #if KB > 217
               movaps 10320(pB), rB2
            #elif KB == 217
               prefA((pfA))
            #endif

         #endif
         #if KB > 217
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10320(pA), rB0
            addps rB0, rC20
            #if KB > 218
               movaps 10336(pB), rB0
            #elif KB == 218
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10320(pA), rB1
            addps rB1, rC21
            #if KB > 218
               movaps 10352(pB), rB1
            #elif KB == 218
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 218
               movaps 10336(pA), rA0
            #elif KB == 218
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 218
               movaps 10352(pA), rA1
            #elif KB == 218
               prefA(-128(pfA))
            #endif

            mulps 10320(pA), rB2
            addps rB2, rC22
            #if KB > 218
               movaps 10368(pB), rB2
            #elif KB == 218
               prefA((pfA))
            #endif

         #endif
         #if KB > 218
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10368(pA), rB0
            addps rB0, rC20
            #if KB > 219
               movaps 10384(pB), rB0
            #elif KB == 219
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10368(pA), rB1
            addps rB1, rC21
            #if KB > 219
               movaps 10400(pB), rB1
            #elif KB == 219
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 219
               movaps 10384(pA), rA0
            #elif KB == 219
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 219
               movaps 10400(pA), rA1
            #elif KB == 219
               prefA(-128(pfA))
            #endif

            mulps 10368(pA), rB2
            addps rB2, rC22
            #if KB > 219
               movaps 10416(pB), rB2
            #elif KB == 219
               prefA((pfA))
            #endif

         #endif
         #if KB > 219
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10416(pA), rB0
            addps rB0, rC20
            #if KB > 220
               movaps 10432(pB), rB0
            #elif KB == 220
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10416(pA), rB1
            addps rB1, rC21
            #if KB > 220
               movaps 10448(pB), rB1
            #elif KB == 220
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 220
               movaps 10432(pA), rA0
            #elif KB == 220
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 220
               movaps 10448(pA), rA1
            #elif KB == 220
               prefA(-128(pfA))
            #endif

            mulps 10416(pA), rB2
            addps rB2, rC22
            #if KB > 220
               movaps 10464(pB), rB2
            #elif KB == 220
               prefA((pfA))
            #endif

         #endif
         #if KB > 220
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10464(pA), rB0
            addps rB0, rC20
            #if KB > 221
               movaps 10480(pB), rB0
            #elif KB == 221
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10464(pA), rB1
            addps rB1, rC21
            #if KB > 221
               movaps 10496(pB), rB1
            #elif KB == 221
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 221
               movaps 10480(pA), rA0
            #elif KB == 221
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 221
               movaps 10496(pA), rA1
            #elif KB == 221
               prefA(-128(pfA))
            #endif

            mulps 10464(pA), rB2
            addps rB2, rC22
            #if KB > 221
               movaps 10512(pB), rB2
            #elif KB == 221
               prefA((pfA))
            #endif

         #endif
         #if KB > 221
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10512(pA), rB0
            addps rB0, rC20
            #if KB > 222
               movaps 10528(pB), rB0
            #elif KB == 222
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10512(pA), rB1
            addps rB1, rC21
            #if KB > 222
               movaps 10544(pB), rB1
            #elif KB == 222
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 222
               movaps 10528(pA), rA0
            #elif KB == 222
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 222
               movaps 10544(pA), rA1
            #elif KB == 222
               prefA(-128(pfA))
            #endif

            mulps 10512(pA), rB2
            addps rB2, rC22
            #if KB > 222
               movaps 10560(pB), rB2
            #elif KB == 222
               prefA((pfA))
            #endif

         #endif
         #if KB > 222
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10560(pA), rB0
            addps rB0, rC20
            #if KB > 223
               movaps 10576(pB), rB0
            #elif KB == 223
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10560(pA), rB1
            addps rB1, rC21
            #if KB > 223
               movaps 10592(pB), rB1
            #elif KB == 223
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 223
               movaps 10576(pA), rA0
            #elif KB == 223
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 223
               movaps 10592(pA), rA1
            #elif KB == 223
               prefA(-128(pfA))
            #endif

            mulps 10560(pA), rB2
            addps rB2, rC22
            #if KB > 223
               movaps 10608(pB), rB2
            #elif KB == 223
               prefA((pfA))
            #endif

         #endif
         #if KB > 223
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10608(pA), rB0
            addps rB0, rC20
            #if KB > 224
               movaps 10624(pB), rB0
            #elif KB == 224
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10608(pA), rB1
            addps rB1, rC21
            #if KB > 224
               movaps 10640(pB), rB1
            #elif KB == 224
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 224
               movaps 10624(pA), rA0
            #elif KB == 224
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 224
               movaps 10640(pA), rA1
            #elif KB == 224
               prefA(-128(pfA))
            #endif

            mulps 10608(pA), rB2
            addps rB2, rC22
            #if KB > 224
               movaps 10656(pB), rB2
            #elif KB == 224
               prefA((pfA))
            #endif

         #endif
         #if KB > 224
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10656(pA), rB0
            addps rB0, rC20
            #if KB > 225
               movaps 10672(pB), rB0
            #elif KB == 225
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10656(pA), rB1
            addps rB1, rC21
            #if KB > 225
               movaps 10688(pB), rB1
            #elif KB == 225
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 225
               movaps 10672(pA), rA0
            #elif KB == 225
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 225
               movaps 10688(pA), rA1
            #elif KB == 225
               prefA(-128(pfA))
            #endif

            mulps 10656(pA), rB2
            addps rB2, rC22
            #if KB > 225
               movaps 10704(pB), rB2
            #elif KB == 225
               prefA((pfA))
            #endif

         #endif
         #if KB > 225
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10704(pA), rB0
            addps rB0, rC20
            #if KB > 226
               movaps 10720(pB), rB0
            #elif KB == 226
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10704(pA), rB1
            addps rB1, rC21
            #if KB > 226
               movaps 10736(pB), rB1
            #elif KB == 226
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 226
               movaps 10720(pA), rA0
            #elif KB == 226
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 226
               movaps 10736(pA), rA1
            #elif KB == 226
               prefA(-128(pfA))
            #endif

            mulps 10704(pA), rB2
            addps rB2, rC22
            #if KB > 226
               movaps 10752(pB), rB2
            #elif KB == 226
               prefA((pfA))
            #endif

         #endif
         #if KB > 226
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10752(pA), rB0
            addps rB0, rC20
            #if KB > 227
               movaps 10768(pB), rB0
            #elif KB == 227
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10752(pA), rB1
            addps rB1, rC21
            #if KB > 227
               movaps 10784(pB), rB1
            #elif KB == 227
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 227
               movaps 10768(pA), rA0
            #elif KB == 227
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 227
               movaps 10784(pA), rA1
            #elif KB == 227
               prefA(-128(pfA))
            #endif

            mulps 10752(pA), rB2
            addps rB2, rC22
            #if KB > 227
               movaps 10800(pB), rB2
            #elif KB == 227
               prefA((pfA))
            #endif

         #endif
         #if KB > 227
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10800(pA), rB0
            addps rB0, rC20
            #if KB > 228
               movaps 10816(pB), rB0
            #elif KB == 228
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10800(pA), rB1
            addps rB1, rC21
            #if KB > 228
               movaps 10832(pB), rB1
            #elif KB == 228
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 228
               movaps 10816(pA), rA0
            #elif KB == 228
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 228
               movaps 10832(pA), rA1
            #elif KB == 228
               prefA(-128(pfA))
            #endif

            mulps 10800(pA), rB2
            addps rB2, rC22
            #if KB > 228
               movaps 10848(pB), rB2
            #elif KB == 228
               prefA((pfA))
            #endif

         #endif
         #if KB > 228
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10848(pA), rB0
            addps rB0, rC20
            #if KB > 229
               movaps 10864(pB), rB0
            #elif KB == 229
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10848(pA), rB1
            addps rB1, rC21
            #if KB > 229
               movaps 10880(pB), rB1
            #elif KB == 229
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 229
               movaps 10864(pA), rA0
            #elif KB == 229
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 229
               movaps 10880(pA), rA1
            #elif KB == 229
               prefA(-128(pfA))
            #endif

            mulps 10848(pA), rB2
            addps rB2, rC22
            #if KB > 229
               movaps 10896(pB), rB2
            #elif KB == 229
               prefA((pfA))
            #endif

         #endif
         #if KB > 229
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10896(pA), rB0
            addps rB0, rC20
            #if KB > 230
               movaps 10912(pB), rB0
            #elif KB == 230
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10896(pA), rB1
            addps rB1, rC21
            #if KB > 230
               movaps 10928(pB), rB1
            #elif KB == 230
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 230
               movaps 10912(pA), rA0
            #elif KB == 230
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 230
               movaps 10928(pA), rA1
            #elif KB == 230
               prefA(-128(pfA))
            #endif

            mulps 10896(pA), rB2
            addps rB2, rC22
            #if KB > 230
               movaps 10944(pB), rB2
            #elif KB == 230
               prefA((pfA))
            #endif

         #endif
         #if KB > 230
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10944(pA), rB0
            addps rB0, rC20
            #if KB > 231
               movaps 10960(pB), rB0
            #elif KB == 231
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10944(pA), rB1
            addps rB1, rC21
            #if KB > 231
               movaps 10976(pB), rB1
            #elif KB == 231
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 231
               movaps 10960(pA), rA0
            #elif KB == 231
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 231
               movaps 10976(pA), rA1
            #elif KB == 231
               prefA(-128(pfA))
            #endif

            mulps 10944(pA), rB2
            addps rB2, rC22
            #if KB > 231
               movaps 10992(pB), rB2
            #elif KB == 231
               prefA((pfA))
            #endif

         #endif
         #if KB > 231
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 10992(pA), rB0
            addps rB0, rC20
            #if KB > 232
               movaps 11008(pB), rB0
            #elif KB == 232
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 10992(pA), rB1
            addps rB1, rC21
            #if KB > 232
               movaps 11024(pB), rB1
            #elif KB == 232
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 232
               movaps 11008(pA), rA0
            #elif KB == 232
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 232
               movaps 11024(pA), rA1
            #elif KB == 232
               prefA(-128(pfA))
            #endif

            mulps 10992(pA), rB2
            addps rB2, rC22
            #if KB > 232
               movaps 11040(pB), rB2
            #elif KB == 232
               prefA((pfA))
            #endif

         #endif
         #if KB > 232
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11040(pA), rB0
            addps rB0, rC20
            #if KB > 233
               movaps 11056(pB), rB0
            #elif KB == 233
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11040(pA), rB1
            addps rB1, rC21
            #if KB > 233
               movaps 11072(pB), rB1
            #elif KB == 233
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 233
               movaps 11056(pA), rA0
            #elif KB == 233
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 233
               movaps 11072(pA), rA1
            #elif KB == 233
               prefA(-128(pfA))
            #endif

            mulps 11040(pA), rB2
            addps rB2, rC22
            #if KB > 233
               movaps 11088(pB), rB2
            #elif KB == 233
               prefA((pfA))
            #endif

         #endif
         #if KB > 233
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11088(pA), rB0
            addps rB0, rC20
            #if KB > 234
               movaps 11104(pB), rB0
            #elif KB == 234
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11088(pA), rB1
            addps rB1, rC21
            #if KB > 234
               movaps 11120(pB), rB1
            #elif KB == 234
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 234
               movaps 11104(pA), rA0
            #elif KB == 234
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 234
               movaps 11120(pA), rA1
            #elif KB == 234
               prefA(-128(pfA))
            #endif

            mulps 11088(pA), rB2
            addps rB2, rC22
            #if KB > 234
               movaps 11136(pB), rB2
            #elif KB == 234
               prefA((pfA))
            #endif

         #endif
         #if KB > 234
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11136(pA), rB0
            addps rB0, rC20
            #if KB > 235
               movaps 11152(pB), rB0
            #elif KB == 235
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11136(pA), rB1
            addps rB1, rC21
            #if KB > 235
               movaps 11168(pB), rB1
            #elif KB == 235
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 235
               movaps 11152(pA), rA0
            #elif KB == 235
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 235
               movaps 11168(pA), rA1
            #elif KB == 235
               prefA(-128(pfA))
            #endif

            mulps 11136(pA), rB2
            addps rB2, rC22
            #if KB > 235
               movaps 11184(pB), rB2
            #elif KB == 235
               prefA((pfA))
            #endif

         #endif
         #if KB > 235
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11184(pA), rB0
            addps rB0, rC20
            #if KB > 236
               movaps 11200(pB), rB0
            #elif KB == 236
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11184(pA), rB1
            addps rB1, rC21
            #if KB > 236
               movaps 11216(pB), rB1
            #elif KB == 236
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 236
               movaps 11200(pA), rA0
            #elif KB == 236
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 236
               movaps 11216(pA), rA1
            #elif KB == 236
               prefA(-128(pfA))
            #endif

            mulps 11184(pA), rB2
            addps rB2, rC22
            #if KB > 236
               movaps 11232(pB), rB2
            #elif KB == 236
               prefA((pfA))
            #endif

         #endif
         #if KB > 236
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11232(pA), rB0
            addps rB0, rC20
            #if KB > 237
               movaps 11248(pB), rB0
            #elif KB == 237
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11232(pA), rB1
            addps rB1, rC21
            #if KB > 237
               movaps 11264(pB), rB1
            #elif KB == 237
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 237
               movaps 11248(pA), rA0
            #elif KB == 237
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 237
               movaps 11264(pA), rA1
            #elif KB == 237
               prefA(-128(pfA))
            #endif

            mulps 11232(pA), rB2
            addps rB2, rC22
            #if KB > 237
               movaps 11280(pB), rB2
            #elif KB == 237
               prefA((pfA))
            #endif

         #endif
         #if KB > 237
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11280(pA), rB0
            addps rB0, rC20
            #if KB > 238
               movaps 11296(pB), rB0
            #elif KB == 238
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11280(pA), rB1
            addps rB1, rC21
            #if KB > 238
               movaps 11312(pB), rB1
            #elif KB == 238
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 238
               movaps 11296(pA), rA0
            #elif KB == 238
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 238
               movaps 11312(pA), rA1
            #elif KB == 238
               prefA(-128(pfA))
            #endif

            mulps 11280(pA), rB2
            addps rB2, rC22
            #if KB > 238
               movaps 11328(pB), rB2
            #elif KB == 238
               prefA((pfA))
            #endif

         #endif
         #if KB > 238
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11328(pA), rB0
            addps rB0, rC20
            #if KB > 239
               movaps 11344(pB), rB0
            #elif KB == 239
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11328(pA), rB1
            addps rB1, rC21
            #if KB > 239
               movaps 11360(pB), rB1
            #elif KB == 239
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 239
               movaps 11344(pA), rA0
            #elif KB == 239
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 239
               movaps 11360(pA), rA1
            #elif KB == 239
               prefA(-128(pfA))
            #endif

            mulps 11328(pA), rB2
            addps rB2, rC22
            #if KB > 239
               movaps 11376(pB), rB2
            #elif KB == 239
               prefA((pfA))
            #endif

         #endif
         #if KB > 239
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11376(pA), rB0
            addps rB0, rC20
            #if KB > 240
               movaps 11392(pB), rB0
            #elif KB == 240
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11376(pA), rB1
            addps rB1, rC21
            #if KB > 240
               movaps 11408(pB), rB1
            #elif KB == 240
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 240
               movaps 11392(pA), rA0
            #elif KB == 240
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 240
               movaps 11408(pA), rA1
            #elif KB == 240
               prefA(-128(pfA))
            #endif

            mulps 11376(pA), rB2
            addps rB2, rC22
            #if KB > 240
               movaps 11424(pB), rB2
            #elif KB == 240
               prefA((pfA))
            #endif

         #endif
         #if KB > 240
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11424(pA), rB0
            addps rB0, rC20
            #if KB > 241
               movaps 11440(pB), rB0
            #elif KB == 241
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11424(pA), rB1
            addps rB1, rC21
            #if KB > 241
               movaps 11456(pB), rB1
            #elif KB == 241
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 241
               movaps 11440(pA), rA0
            #elif KB == 241
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 241
               movaps 11456(pA), rA1
            #elif KB == 241
               prefA(-128(pfA))
            #endif

            mulps 11424(pA), rB2
            addps rB2, rC22
            #if KB > 241
               movaps 11472(pB), rB2
            #elif KB == 241
               prefA((pfA))
            #endif

         #endif
         #if KB > 241
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11472(pA), rB0
            addps rB0, rC20
            #if KB > 242
               movaps 11488(pB), rB0
            #elif KB == 242
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11472(pA), rB1
            addps rB1, rC21
            #if KB > 242
               movaps 11504(pB), rB1
            #elif KB == 242
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 242
               movaps 11488(pA), rA0
            #elif KB == 242
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 242
               movaps 11504(pA), rA1
            #elif KB == 242
               prefA(-128(pfA))
            #endif

            mulps 11472(pA), rB2
            addps rB2, rC22
            #if KB > 242
               movaps 11520(pB), rB2
            #elif KB == 242
               prefA((pfA))
            #endif

         #endif
         #if KB > 242
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11520(pA), rB0
            addps rB0, rC20
            #if KB > 243
               movaps 11536(pB), rB0
            #elif KB == 243
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11520(pA), rB1
            addps rB1, rC21
            #if KB > 243
               movaps 11552(pB), rB1
            #elif KB == 243
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 243
               movaps 11536(pA), rA0
            #elif KB == 243
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 243
               movaps 11552(pA), rA1
            #elif KB == 243
               prefA(-128(pfA))
            #endif

            mulps 11520(pA), rB2
            addps rB2, rC22
            #if KB > 243
               movaps 11568(pB), rB2
            #elif KB == 243
               prefA((pfA))
            #endif

         #endif
         #if KB > 243
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11568(pA), rB0
            addps rB0, rC20
            #if KB > 244
               movaps 11584(pB), rB0
            #elif KB == 244
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11568(pA), rB1
            addps rB1, rC21
            #if KB > 244
               movaps 11600(pB), rB1
            #elif KB == 244
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 244
               movaps 11584(pA), rA0
            #elif KB == 244
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 244
               movaps 11600(pA), rA1
            #elif KB == 244
               prefA(-128(pfA))
            #endif

            mulps 11568(pA), rB2
            addps rB2, rC22
            #if KB > 244
               movaps 11616(pB), rB2
            #elif KB == 244
               prefA((pfA))
            #endif

         #endif
         #if KB > 244
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11616(pA), rB0
            addps rB0, rC20
            #if KB > 245
               movaps 11632(pB), rB0
            #elif KB == 245
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11616(pA), rB1
            addps rB1, rC21
            #if KB > 245
               movaps 11648(pB), rB1
            #elif KB == 245
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 245
               movaps 11632(pA), rA0
            #elif KB == 245
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 245
               movaps 11648(pA), rA1
            #elif KB == 245
               prefA(-128(pfA))
            #endif

            mulps 11616(pA), rB2
            addps rB2, rC22
            #if KB > 245
               movaps 11664(pB), rB2
            #elif KB == 245
               prefA((pfA))
            #endif

         #endif
         #if KB > 245
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11664(pA), rB0
            addps rB0, rC20
            #if KB > 246
               movaps 11680(pB), rB0
            #elif KB == 246
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11664(pA), rB1
            addps rB1, rC21
            #if KB > 246
               movaps 11696(pB), rB1
            #elif KB == 246
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 246
               movaps 11680(pA), rA0
            #elif KB == 246
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 246
               movaps 11696(pA), rA1
            #elif KB == 246
               prefA(-128(pfA))
            #endif

            mulps 11664(pA), rB2
            addps rB2, rC22
            #if KB > 246
               movaps 11712(pB), rB2
            #elif KB == 246
               prefA((pfA))
            #endif

         #endif
         #if KB > 246
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11712(pA), rB0
            addps rB0, rC20
            #if KB > 247
               movaps 11728(pB), rB0
            #elif KB == 247
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11712(pA), rB1
            addps rB1, rC21
            #if KB > 247
               movaps 11744(pB), rB1
            #elif KB == 247
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 247
               movaps 11728(pA), rA0
            #elif KB == 247
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 247
               movaps 11744(pA), rA1
            #elif KB == 247
               prefA(-128(pfA))
            #endif

            mulps 11712(pA), rB2
            addps rB2, rC22
            #if KB > 247
               movaps 11760(pB), rB2
            #elif KB == 247
               prefA((pfA))
            #endif

         #endif
         #if KB > 247
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11760(pA), rB0
            addps rB0, rC20
            #if KB > 248
               movaps 11776(pB), rB0
            #elif KB == 248
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11760(pA), rB1
            addps rB1, rC21
            #if KB > 248
               movaps 11792(pB), rB1
            #elif KB == 248
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 248
               movaps 11776(pA), rA0
            #elif KB == 248
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 248
               movaps 11792(pA), rA1
            #elif KB == 248
               prefA(-128(pfA))
            #endif

            mulps 11760(pA), rB2
            addps rB2, rC22
            #if KB > 248
               movaps 11808(pB), rB2
            #elif KB == 248
               prefA((pfA))
            #endif

         #endif
         #if KB > 248
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11808(pA), rB0
            addps rB0, rC20
            #if KB > 249
               movaps 11824(pB), rB0
            #elif KB == 249
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11808(pA), rB1
            addps rB1, rC21
            #if KB > 249
               movaps 11840(pB), rB1
            #elif KB == 249
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 249
               movaps 11824(pA), rA0
            #elif KB == 249
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 249
               movaps 11840(pA), rA1
            #elif KB == 249
               prefA(-128(pfA))
            #endif

            mulps 11808(pA), rB2
            addps rB2, rC22
            #if KB > 249
               movaps 11856(pB), rB2
            #elif KB == 249
               prefA((pfA))
            #endif

         #endif
         #if KB > 249
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11856(pA), rB0
            addps rB0, rC20
            #if KB > 250
               movaps 11872(pB), rB0
            #elif KB == 250
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11856(pA), rB1
            addps rB1, rC21
            #if KB > 250
               movaps 11888(pB), rB1
            #elif KB == 250
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 250
               movaps 11872(pA), rA0
            #elif KB == 250
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 250
               movaps 11888(pA), rA1
            #elif KB == 250
               prefA(-128(pfA))
            #endif

            mulps 11856(pA), rB2
            addps rB2, rC22
            #if KB > 250
               movaps 11904(pB), rB2
            #elif KB == 250
               prefA((pfA))
            #endif

         #endif
         #if KB > 250
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11904(pA), rB0
            addps rB0, rC20
            #if KB > 251
               movaps 11920(pB), rB0
            #elif KB == 251
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11904(pA), rB1
            addps rB1, rC21
            #if KB > 251
               movaps 11936(pB), rB1
            #elif KB == 251
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 251
               movaps 11920(pA), rA0
            #elif KB == 251
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 251
               movaps 11936(pA), rA1
            #elif KB == 251
               prefA(-128(pfA))
            #endif

            mulps 11904(pA), rB2
            addps rB2, rC22
            #if KB > 251
               movaps 11952(pB), rB2
            #elif KB == 251
               prefA((pfA))
            #endif

         #endif
         #if KB > 251
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 11952(pA), rB0
            addps rB0, rC20
            #if KB > 252
               movaps 11968(pB), rB0
            #elif KB == 252
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 11952(pA), rB1
            addps rB1, rC21
            #if KB > 252
               movaps 11984(pB), rB1
            #elif KB == 252
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 252
               movaps 11968(pA), rA0
            #elif KB == 252
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 252
               movaps 11984(pA), rA1
            #elif KB == 252
               prefA(-128(pfA))
            #endif

            mulps 11952(pA), rB2
            addps rB2, rC22
            #if KB > 252
               movaps 12000(pB), rB2
            #elif KB == 252
               prefA((pfA))
            #endif

         #endif
         #if KB > 252
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 12000(pA), rB0
            addps rB0, rC20
            #if KB > 253
               movaps 12016(pB), rB0
            #elif KB == 253
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 12000(pA), rB1
            addps rB1, rC21
            #if KB > 253
               movaps 12032(pB), rB1
            #elif KB == 253
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 253
               movaps 12016(pA), rA0
            #elif KB == 253
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 253
               movaps 12032(pA), rA1
            #elif KB == 253
               prefA(-128(pfA))
            #endif

            mulps 12000(pA), rB2
            addps rB2, rC22
            #if KB > 253
               movaps 12048(pB), rB2
            #elif KB == 253
               prefA((pfA))
            #endif

         #endif
         #if KB > 253
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 12048(pA), rB0
            addps rB0, rC20
            #if KB > 254
               movaps 12064(pB), rB0
            #elif KB == 254
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 12048(pA), rB1
            addps rB1, rC21
            #if KB > 254
               movaps 12080(pB), rB1
            #elif KB == 254
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 254
               movaps 12064(pA), rA0
            #elif KB == 254
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 254
               movaps 12080(pA), rA1
            #elif KB == 254
               prefA(-128(pfA))
            #endif

            mulps 12048(pA), rB2
            addps rB2, rC22
            #if KB > 254
               movaps 12096(pB), rB2
            #elif KB == 254
               prefA((pfA))
            #endif

         #endif
         #if KB > 254
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 12096(pA), rB0
            addps rB0, rC20
            #if KB > 255
               movaps 12112(pB), rB0
            #elif KB == 255
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 12096(pA), rB1
            addps rB1, rC21
            #if KB > 255
               movaps 12128(pB), rB1
            #elif KB == 255
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 255
               movaps 12112(pA), rA0
            #elif KB == 255
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 255
               movaps 12128(pA), rA1
            #elif KB == 255
               prefA(-128(pfA))
            #endif

            mulps 12096(pA), rB2
            addps rB2, rC22
            #if KB > 255
               movaps 12144(pB), rB2
            #elif KB == 255
               prefA((pfA))
            #endif

         #endif
         #if KB > 255
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            mulps 12144(pA), rB0
            addps rB0, rC20
            #if KB > 256
               movaps 12160(pB), rB0
            #elif KB == 256
               prefB(-128(pfB))
            #endif

            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01

            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11

            mulps 12144(pA), rB1
            addps rB1, rC21
            #if KB > 256
               movaps 12176(pB), rB1
            #elif KB == 256
               prefB(-64(pfB))
            #endif

            mulps rB2, rA0
            addps rA0, rC02
            #if KB > 256
               movaps 12160(pA), rA0
            #elif KB == 256
               prefB((pfB))
            #endif

            mulps rB2, rA1
            addps rA1, rC12
            #if KB > 256
               movaps 12176(pA), rA1
            #elif KB == 256
               prefA(-128(pfA))
            #endif

            mulps 12144(pA), rB2
            addps rB2, rC22
            #if KB > 256
               movaps 12192(pB), rB2
            #elif KB == 256
               prefA((pfA))
            #endif

         #endif

         add incPF, pfA
         add incPF, pfB
         #ifndef BETA0
            VCOP -128(pC), rC00
         #endif
         movaps rC00, -128(pC)
         #ifndef BETA0
            VCOP -112(pC), rC10
         #endif
         movaps rC10, -112(pC)
         #ifndef BETA0
            VCOP -96(pC), rC20
         #endif
         movaps rC20, -96(pC)
         #ifndef BETA0
            VCOP -80(pC), rC01
         #endif
         movaps rC01, -80(pC)
         #ifndef BETA0
            VCOP -64(pC), rC11
         #endif
         movaps rC11, -64(pC)
         #ifndef BETA0
            VCOP -48(pC), rC21
         #endif
         movaps rC21, -48(pC)
         #ifndef BETA0
            VCOP -32(pC), rC02
         #endif
         movaps rC02, -32(pC)
         #ifndef BETA0
            VCOP -16(pC), rC12
         #endif
         movaps rC12, -16(pC)
         #ifndef BETA0
            VCOP (pC), rC22
         #endif
         movaps rC22, (pC)
         add $KB*3*16, pB
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
