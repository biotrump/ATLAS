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
#if KB < 45
   #define prefA(m_)
   #define prefB(m_) prefetcht0 m_
#else
   #if KB > 78
      #define prefA(m_)
   #else
      #define prefA(m_) prefetcht2 m_
   #endif
   #define prefB(m_) prefetcht2 m_
#endif
#ifdef ATL_3DNow
   #define prefC(m_) prefetchw m_
#else
   #define prefC(m_) prefetcht0 m_
#endif
#if defined(BETAN) || defined(BETAn)
   #define BETAN1
#endif
#ifdef BETAN1
   #define VCOP subpd
#else
   #define VCOP addpd
#endif
/* #define movapd movaps */
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
   #if KB > 78  /* only pfB is prefetched, so put most imp mat in it */
      #if defined(ATL_MOVEB)
         lea KB*6*8(pB), pfB
      #else
         lea KB*12*8(pA), pfB
      #endif
   #else
      movq 8(%rsp), pfA       /* pfA = pAn */
      movq 16(%rsp), pfB      /* pfB = pBn */
      cmp pfB, pB
      CMOVE pfA, pfB
      CMOVE pA, pfA
   #endif
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
         movddup -128(pB), rC00
         movapd -128(pA), rA0
         movapd rC00, rC10

         mulpd rA0, rC00
         movapd rC10, rC20
         movapd -112(pA), rA1

         mulpd rA1, rC10
         movddup -120(pB), rC01
         movapd -96(pA), rA2

         mulpd rA2, rC20
         movapd rC01, rC11
         movddup -112(pB), rC02

         mulpd rA0, rC01
         movapd rC11, rC21
         movddup -104(pB), rB0

         movapd rC02, rC12
         mulpd rA1, rC11
         movddup -96(pB), rB1

         mulpd rA2, rC21
         movapd rC12, rC22
         prefC(-128(pC))

         mulpd rA0, rC02
         movapd -80(pA), rA0
         prefC(-64(pC))

         mulpd rA1, rC12
         movapd -64(pA), rA1
         movddup -88(pB), rB2

         mulpd rA2, rC22
         prefC((pC))
         prefC(64(pC))

         #if KB > 1
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00
            #if defined(ATL_MOVEC)
               prefC(128-192(pC,r24))
            #endif

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10
            #if defined(ATL_MOVEC)
               prefC((pC,r24))
            #endif

            mulpd -48(pA), rB0
            addpd rB0, rC20
            #if KB > 2
               movddup -80(pB), rB0
            #elif KB == 2
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd -48(pA), rB1
            addpd rB1, rC21
            #if KB > 2
               movddup -72(pB), rB1
            #elif KB == 2
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 2
               movaps -32(pA), rA0
            #elif KB == 2
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 2
               movaps -16(pA), rA1
            #elif KB == 2
               prefA(-128(pfA))
            #endif

            mulpd -48(pA), rB2
            addpd rB2, rC22
            #if KB > 2
               movddup -64(pB), rB2
            #elif KB == 2
               prefA((pfA))
            #endif

         #endif
         #if KB > 2
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 0(pA), rB0
            addpd rB0, rC20
            #if KB > 3
               movddup -56(pB), rB0
            #elif KB == 3
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 0(pA), rB1
            addpd rB1, rC21
            #if KB > 3
               movddup -48(pB), rB1
            #elif KB == 3
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 3
               movaps 16(pA), rA0
            #elif KB == 3
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 3
               movaps 32(pA), rA1
            #elif KB == 3
               prefA(-128(pfA))
            #endif

            mulpd 0(pA), rB2
            addpd rB2, rC22
            #if KB > 3
               movddup -40(pB), rB2
            #elif KB == 3
               prefA((pfA))
            #endif

         #endif
         #if KB > 3
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 48(pA), rB0
            addpd rB0, rC20
            #if KB > 4
               movddup -32(pB), rB0
            #elif KB == 4
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 48(pA), rB1
            addpd rB1, rC21
            #if KB > 4
               movddup -24(pB), rB1
            #elif KB == 4
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 4
               movaps 64(pA), rA0
            #elif KB == 4
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 4
               movaps 80(pA), rA1
            #elif KB == 4
               prefA(-128(pfA))
            #endif

            mulpd 48(pA), rB2
            addpd rB2, rC22
            #if KB > 4
               movddup -16(pB), rB2
            #elif KB == 4
               prefA((pfA))
            #endif

         #endif
         #if KB > 4
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 96(pA), rB0
            addpd rB0, rC20
            #if KB > 5
               movddup -8(pB), rB0
            #elif KB == 5
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 96(pA), rB1
            addpd rB1, rC21
            #if KB > 5
               movddup 0(pB), rB1
            #elif KB == 5
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 5
               movaps 112(pA), rA0
            #elif KB == 5
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 5
               movaps 128(pA), rA1
            #elif KB == 5
               prefA(-128(pfA))
            #endif

            mulpd 96(pA), rB2
            addpd rB2, rC22
            #if KB > 5
               movddup 8(pB), rB2
            #elif KB == 5
               prefA((pfA))
            #endif

         #endif
         #if KB > 5
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 144(pA), rB0
            addpd rB0, rC20
            #if KB > 6
               movddup 16(pB), rB0
            #elif KB == 6
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 144(pA), rB1
            addpd rB1, rC21
            #if KB > 6
               movddup 24(pB), rB1
            #elif KB == 6
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 6
               movaps 160(pA), rA0
            #elif KB == 6
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 6
               movaps 176(pA), rA1
            #elif KB == 6
               prefA(-128(pfA))
            #endif

            mulpd 144(pA), rB2
            addpd rB2, rC22
            #if KB > 6
               movddup 32(pB), rB2
            #elif KB == 6
               prefA((pfA))
            #endif

         #endif
         #if KB > 6
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 192(pA), rB0
            addpd rB0, rC20
            #if KB > 7
               movddup 40(pB), rB0
            #elif KB == 7
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 192(pA), rB1
            addpd rB1, rC21
            #if KB > 7
               movddup 48(pB), rB1
            #elif KB == 7
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 7
               movaps 208(pA), rA0
            #elif KB == 7
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 7
               movaps 224(pA), rA1
            #elif KB == 7
               prefA(-128(pfA))
            #endif

            mulpd 192(pA), rB2
            addpd rB2, rC22
            #if KB > 7
               movddup 56(pB), rB2
            #elif KB == 7
               prefA((pfA))
            #endif

         #endif
         #if KB > 7
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 240(pA), rB0
            addpd rB0, rC20
            #if KB > 8
               movddup 64(pB), rB0
            #elif KB == 8
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 240(pA), rB1
            addpd rB1, rC21
            #if KB > 8
               movddup 72(pB), rB1
            #elif KB == 8
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 8
               movaps 256(pA), rA0
            #elif KB == 8
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 8
               movaps 272(pA), rA1
            #elif KB == 8
               prefA(-128(pfA))
            #endif

            mulpd 240(pA), rB2
            addpd rB2, rC22
            #if KB > 8
               movddup 80(pB), rB2
            #elif KB == 8
               prefA((pfA))
            #endif

         #endif
         #if KB > 8
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 288(pA), rB0
            addpd rB0, rC20
            #if KB > 9
               movddup 88(pB), rB0
            #elif KB == 9
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 288(pA), rB1
            addpd rB1, rC21
            #if KB > 9
               movddup 96(pB), rB1
            #elif KB == 9
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 9
               movaps 304(pA), rA0
            #elif KB == 9
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 9
               movaps 320(pA), rA1
            #elif KB == 9
               prefA(-128(pfA))
            #endif

            mulpd 288(pA), rB2
            addpd rB2, rC22
            #if KB > 9
               movddup 104(pB), rB2
            #elif KB == 9
               prefA((pfA))
            #endif

         #endif
         #if KB > 9
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 336(pA), rB0
            addpd rB0, rC20
            #if KB > 10
               movddup 112(pB), rB0
            #elif KB == 10
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 336(pA), rB1
            addpd rB1, rC21
            #if KB > 10
               movddup 120(pB), rB1
            #elif KB == 10
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 10
               movaps 352(pA), rA0
            #elif KB == 10
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 10
               movaps 368(pA), rA1
            #elif KB == 10
               prefA(-128(pfA))
            #endif

            mulpd 336(pA), rB2
            addpd rB2, rC22
            #if KB > 10
               movddup 128(pB), rB2
            #elif KB == 10
               prefA((pfA))
            #endif

         #endif
         #if KB > 10
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 384(pA), rB0
            addpd rB0, rC20
            #if KB > 11
               movddup 136(pB), rB0
            #elif KB == 11
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 384(pA), rB1
            addpd rB1, rC21
            #if KB > 11
               movddup 144(pB), rB1
            #elif KB == 11
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 11
               movaps 400(pA), rA0
            #elif KB == 11
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 11
               movaps 416(pA), rA1
            #elif KB == 11
               prefA(-128(pfA))
            #endif

            mulpd 384(pA), rB2
            addpd rB2, rC22
            #if KB > 11
               movddup 152(pB), rB2
            #elif KB == 11
               prefA((pfA))
            #endif

         #endif
         #if KB > 11
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 432(pA), rB0
            addpd rB0, rC20
            #if KB > 12
               movddup 160(pB), rB0
            #elif KB == 12
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 432(pA), rB1
            addpd rB1, rC21
            #if KB > 12
               movddup 168(pB), rB1
            #elif KB == 12
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 12
               movaps 448(pA), rA0
            #elif KB == 12
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 12
               movaps 464(pA), rA1
            #elif KB == 12
               prefA(-128(pfA))
            #endif

            mulpd 432(pA), rB2
            addpd rB2, rC22
            #if KB > 12
               movddup 176(pB), rB2
            #elif KB == 12
               prefA((pfA))
            #endif

         #endif
         #if KB > 12
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 480(pA), rB0
            addpd rB0, rC20
            #if KB > 13
               movddup 184(pB), rB0
            #elif KB == 13
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 480(pA), rB1
            addpd rB1, rC21
            #if KB > 13
               movddup 192(pB), rB1
            #elif KB == 13
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 13
               movaps 496(pA), rA0
            #elif KB == 13
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 13
               movaps 512(pA), rA1
            #elif KB == 13
               prefA(-128(pfA))
            #endif

            mulpd 480(pA), rB2
            addpd rB2, rC22
            #if KB > 13
               movddup 200(pB), rB2
            #elif KB == 13
               prefA((pfA))
            #endif

         #endif
         #if KB > 13
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 528(pA), rB0
            addpd rB0, rC20
            #if KB > 14
               movddup 208(pB), rB0
            #elif KB == 14
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 528(pA), rB1
            addpd rB1, rC21
            #if KB > 14
               movddup 216(pB), rB1
            #elif KB == 14
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 14
               movaps 544(pA), rA0
            #elif KB == 14
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 14
               movaps 560(pA), rA1
            #elif KB == 14
               prefA(-128(pfA))
            #endif

            mulpd 528(pA), rB2
            addpd rB2, rC22
            #if KB > 14
               movddup 224(pB), rB2
            #elif KB == 14
               prefA((pfA))
            #endif

         #endif
         #if KB > 14
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 576(pA), rB0
            addpd rB0, rC20
            #if KB > 15
               movddup 232(pB), rB0
            #elif KB == 15
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 576(pA), rB1
            addpd rB1, rC21
            #if KB > 15
               movddup 240(pB), rB1
            #elif KB == 15
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 15
               movaps 592(pA), rA0
            #elif KB == 15
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 15
               movaps 608(pA), rA1
            #elif KB == 15
               prefA(-128(pfA))
            #endif

            mulpd 576(pA), rB2
            addpd rB2, rC22
            #if KB > 15
               movddup 248(pB), rB2
            #elif KB == 15
               prefA((pfA))
            #endif

         #endif
         #if KB > 15
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 624(pA), rB0
            addpd rB0, rC20
            #if KB > 16
               movddup 256(pB), rB0
            #elif KB == 16
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 624(pA), rB1
            addpd rB1, rC21
            #if KB > 16
               movddup 264(pB), rB1
            #elif KB == 16
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 16
               movaps 640(pA), rA0
            #elif KB == 16
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 16
               movaps 656(pA), rA1
            #elif KB == 16
               prefA(-128(pfA))
            #endif

            mulpd 624(pA), rB2
            addpd rB2, rC22
            #if KB > 16
               movddup 272(pB), rB2
            #elif KB == 16
               prefA((pfA))
            #endif

         #endif
         #if KB > 16
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 672(pA), rB0
            addpd rB0, rC20
            #if KB > 17
               movddup 280(pB), rB0
            #elif KB == 17
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 672(pA), rB1
            addpd rB1, rC21
            #if KB > 17
               movddup 288(pB), rB1
            #elif KB == 17
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 17
               movaps 688(pA), rA0
            #elif KB == 17
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 17
               movaps 704(pA), rA1
            #elif KB == 17
               prefA(-128(pfA))
            #endif

            mulpd 672(pA), rB2
            addpd rB2, rC22
            #if KB > 17
               movddup 296(pB), rB2
            #elif KB == 17
               prefA((pfA))
            #endif

         #endif
         #if KB > 17
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 720(pA), rB0
            addpd rB0, rC20
            #if KB > 18
               movddup 304(pB), rB0
            #elif KB == 18
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 720(pA), rB1
            addpd rB1, rC21
            #if KB > 18
               movddup 312(pB), rB1
            #elif KB == 18
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 18
               movaps 736(pA), rA0
            #elif KB == 18
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 18
               movaps 752(pA), rA1
            #elif KB == 18
               prefA(-128(pfA))
            #endif

            mulpd 720(pA), rB2
            addpd rB2, rC22
            #if KB > 18
               movddup 320(pB), rB2
            #elif KB == 18
               prefA((pfA))
            #endif

         #endif
         #if KB > 18
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 768(pA), rB0
            addpd rB0, rC20
            #if KB > 19
               movddup 328(pB), rB0
            #elif KB == 19
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 768(pA), rB1
            addpd rB1, rC21
            #if KB > 19
               movddup 336(pB), rB1
            #elif KB == 19
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 19
               movaps 784(pA), rA0
            #elif KB == 19
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 19
               movaps 800(pA), rA1
            #elif KB == 19
               prefA(-128(pfA))
            #endif

            mulpd 768(pA), rB2
            addpd rB2, rC22
            #if KB > 19
               movddup 344(pB), rB2
            #elif KB == 19
               prefA((pfA))
            #endif

         #endif
         #if KB > 19
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 816(pA), rB0
            addpd rB0, rC20
            #if KB > 20
               movddup 352(pB), rB0
            #elif KB == 20
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 816(pA), rB1
            addpd rB1, rC21
            #if KB > 20
               movddup 360(pB), rB1
            #elif KB == 20
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 20
               movaps 832(pA), rA0
            #elif KB == 20
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 20
               movaps 848(pA), rA1
            #elif KB == 20
               prefA(-128(pfA))
            #endif

            mulpd 816(pA), rB2
            addpd rB2, rC22
            #if KB > 20
               movddup 368(pB), rB2
            #elif KB == 20
               prefA((pfA))
            #endif

         #endif
         #if KB > 20
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 864(pA), rB0
            addpd rB0, rC20
            #if KB > 21
               movddup 376(pB), rB0
            #elif KB == 21
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 864(pA), rB1
            addpd rB1, rC21
            #if KB > 21
               movddup 384(pB), rB1
            #elif KB == 21
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 21
               movaps 880(pA), rA0
            #elif KB == 21
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 21
               movaps 896(pA), rA1
            #elif KB == 21
               prefA(-128(pfA))
            #endif

            mulpd 864(pA), rB2
            addpd rB2, rC22
            #if KB > 21
               movddup 392(pB), rB2
            #elif KB == 21
               prefA((pfA))
            #endif

         #endif
         #if KB > 21
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 912(pA), rB0
            addpd rB0, rC20
            #if KB > 22
               movddup 400(pB), rB0
            #elif KB == 22
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 912(pA), rB1
            addpd rB1, rC21
            #if KB > 22
               movddup 408(pB), rB1
            #elif KB == 22
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 22
               movaps 928(pA), rA0
            #elif KB == 22
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 22
               movaps 944(pA), rA1
            #elif KB == 22
               prefA(-128(pfA))
            #endif

            mulpd 912(pA), rB2
            addpd rB2, rC22
            #if KB > 22
               movddup 416(pB), rB2
            #elif KB == 22
               prefA((pfA))
            #endif

         #endif
         #if KB > 22
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 960(pA), rB0
            addpd rB0, rC20
            #if KB > 23
               movddup 424(pB), rB0
            #elif KB == 23
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 960(pA), rB1
            addpd rB1, rC21
            #if KB > 23
               movddup 432(pB), rB1
            #elif KB == 23
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 23
               movaps 976(pA), rA0
            #elif KB == 23
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 23
               movaps 992(pA), rA1
            #elif KB == 23
               prefA(-128(pfA))
            #endif

            mulpd 960(pA), rB2
            addpd rB2, rC22
            #if KB > 23
               movddup 440(pB), rB2
            #elif KB == 23
               prefA((pfA))
            #endif

         #endif
         #if KB > 23
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1008(pA), rB0
            addpd rB0, rC20
            #if KB > 24
               movddup 448(pB), rB0
            #elif KB == 24
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1008(pA), rB1
            addpd rB1, rC21
            #if KB > 24
               movddup 456(pB), rB1
            #elif KB == 24
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 24
               movaps 1024(pA), rA0
            #elif KB == 24
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 24
               movaps 1040(pA), rA1
            #elif KB == 24
               prefA(-128(pfA))
            #endif

            mulpd 1008(pA), rB2
            addpd rB2, rC22
            #if KB > 24
               movddup 464(pB), rB2
            #elif KB == 24
               prefA((pfA))
            #endif

         #endif
         #if KB > 24
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1056(pA), rB0
            addpd rB0, rC20
            #if KB > 25
               movddup 472(pB), rB0
            #elif KB == 25
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1056(pA), rB1
            addpd rB1, rC21
            #if KB > 25
               movddup 480(pB), rB1
            #elif KB == 25
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 25
               movaps 1072(pA), rA0
            #elif KB == 25
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 25
               movaps 1088(pA), rA1
            #elif KB == 25
               prefA(-128(pfA))
            #endif

            mulpd 1056(pA), rB2
            addpd rB2, rC22
            #if KB > 25
               movddup 488(pB), rB2
            #elif KB == 25
               prefA((pfA))
            #endif

         #endif
         #if KB > 25
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1104(pA), rB0
            addpd rB0, rC20
            #if KB > 26
               movddup 496(pB), rB0
            #elif KB == 26
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1104(pA), rB1
            addpd rB1, rC21
            #if KB > 26
               movddup 504(pB), rB1
            #elif KB == 26
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 26
               movaps 1120(pA), rA0
            #elif KB == 26
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 26
               movaps 1136(pA), rA1
            #elif KB == 26
               prefA(-128(pfA))
            #endif

            mulpd 1104(pA), rB2
            addpd rB2, rC22
            #if KB > 26
               movddup 512(pB), rB2
            #elif KB == 26
               prefA((pfA))
            #endif

         #endif
         #if KB > 26
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1152(pA), rB0
            addpd rB0, rC20
            #if KB > 27
               movddup 520(pB), rB0
            #elif KB == 27
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1152(pA), rB1
            addpd rB1, rC21
            #if KB > 27
               movddup 528(pB), rB1
            #elif KB == 27
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 27
               movaps 1168(pA), rA0
            #elif KB == 27
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 27
               movaps 1184(pA), rA1
            #elif KB == 27
               prefA(-128(pfA))
            #endif

            mulpd 1152(pA), rB2
            addpd rB2, rC22
            #if KB > 27
               movddup 536(pB), rB2
            #elif KB == 27
               prefA((pfA))
            #endif

         #endif
         #if KB > 27
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1200(pA), rB0
            addpd rB0, rC20
            #if KB > 28
               movddup 544(pB), rB0
            #elif KB == 28
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1200(pA), rB1
            addpd rB1, rC21
            #if KB > 28
               movddup 552(pB), rB1
            #elif KB == 28
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 28
               movaps 1216(pA), rA0
            #elif KB == 28
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 28
               movaps 1232(pA), rA1
            #elif KB == 28
               prefA(-128(pfA))
            #endif

            mulpd 1200(pA), rB2
            addpd rB2, rC22
            #if KB > 28
               movddup 560(pB), rB2
            #elif KB == 28
               prefA((pfA))
            #endif

         #endif
         #if KB > 28
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1248(pA), rB0
            addpd rB0, rC20
            #if KB > 29
               movddup 568(pB), rB0
            #elif KB == 29
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1248(pA), rB1
            addpd rB1, rC21
            #if KB > 29
               movddup 576(pB), rB1
            #elif KB == 29
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 29
               movaps 1264(pA), rA0
            #elif KB == 29
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 29
               movaps 1280(pA), rA1
            #elif KB == 29
               prefA(-128(pfA))
            #endif

            mulpd 1248(pA), rB2
            addpd rB2, rC22
            #if KB > 29
               movddup 584(pB), rB2
            #elif KB == 29
               prefA((pfA))
            #endif

         #endif
         #if KB > 29
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1296(pA), rB0
            addpd rB0, rC20
            #if KB > 30
               movddup 592(pB), rB0
            #elif KB == 30
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1296(pA), rB1
            addpd rB1, rC21
            #if KB > 30
               movddup 600(pB), rB1
            #elif KB == 30
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 30
               movaps 1312(pA), rA0
            #elif KB == 30
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 30
               movaps 1328(pA), rA1
            #elif KB == 30
               prefA(-128(pfA))
            #endif

            mulpd 1296(pA), rB2
            addpd rB2, rC22
            #if KB > 30
               movddup 608(pB), rB2
            #elif KB == 30
               prefA((pfA))
            #endif

         #endif
         #if KB > 30
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1344(pA), rB0
            addpd rB0, rC20
            #if KB > 31
               movddup 616(pB), rB0
            #elif KB == 31
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1344(pA), rB1
            addpd rB1, rC21
            #if KB > 31
               movddup 624(pB), rB1
            #elif KB == 31
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 31
               movaps 1360(pA), rA0
            #elif KB == 31
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 31
               movaps 1376(pA), rA1
            #elif KB == 31
               prefA(-128(pfA))
            #endif

            mulpd 1344(pA), rB2
            addpd rB2, rC22
            #if KB > 31
               movddup 632(pB), rB2
            #elif KB == 31
               prefA((pfA))
            #endif

         #endif
         #if KB > 31
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1392(pA), rB0
            addpd rB0, rC20
            #if KB > 32
               movddup 640(pB), rB0
            #elif KB == 32
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1392(pA), rB1
            addpd rB1, rC21
            #if KB > 32
               movddup 648(pB), rB1
            #elif KB == 32
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 32
               movaps 1408(pA), rA0
            #elif KB == 32
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 32
               movaps 1424(pA), rA1
            #elif KB == 32
               prefA(-128(pfA))
            #endif

            mulpd 1392(pA), rB2
            addpd rB2, rC22
            #if KB > 32
               movddup 656(pB), rB2
            #elif KB == 32
               prefA((pfA))
            #endif

         #endif
         #if KB > 32
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1440(pA), rB0
            addpd rB0, rC20
            #if KB > 33
               movddup 664(pB), rB0
            #elif KB == 33
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1440(pA), rB1
            addpd rB1, rC21
            #if KB > 33
               movddup 672(pB), rB1
            #elif KB == 33
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 33
               movaps 1456(pA), rA0
            #elif KB == 33
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 33
               movaps 1472(pA), rA1
            #elif KB == 33
               prefA(-128(pfA))
            #endif

            mulpd 1440(pA), rB2
            addpd rB2, rC22
            #if KB > 33
               movddup 680(pB), rB2
            #elif KB == 33
               prefA((pfA))
            #endif

         #endif
         #if KB > 33
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1488(pA), rB0
            addpd rB0, rC20
            #if KB > 34
               movddup 688(pB), rB0
            #elif KB == 34
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1488(pA), rB1
            addpd rB1, rC21
            #if KB > 34
               movddup 696(pB), rB1
            #elif KB == 34
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 34
               movaps 1504(pA), rA0
            #elif KB == 34
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 34
               movaps 1520(pA), rA1
            #elif KB == 34
               prefA(-128(pfA))
            #endif

            mulpd 1488(pA), rB2
            addpd rB2, rC22
            #if KB > 34
               movddup 704(pB), rB2
            #elif KB == 34
               prefA((pfA))
            #endif

         #endif
         #if KB > 34
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1536(pA), rB0
            addpd rB0, rC20
            #if KB > 35
               movddup 712(pB), rB0
            #elif KB == 35
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1536(pA), rB1
            addpd rB1, rC21
            #if KB > 35
               movddup 720(pB), rB1
            #elif KB == 35
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 35
               movaps 1552(pA), rA0
            #elif KB == 35
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 35
               movaps 1568(pA), rA1
            #elif KB == 35
               prefA(-128(pfA))
            #endif

            mulpd 1536(pA), rB2
            addpd rB2, rC22
            #if KB > 35
               movddup 728(pB), rB2
            #elif KB == 35
               prefA((pfA))
            #endif

         #endif
         #if KB > 35
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1584(pA), rB0
            addpd rB0, rC20
            #if KB > 36
               movddup 736(pB), rB0
            #elif KB == 36
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1584(pA), rB1
            addpd rB1, rC21
            #if KB > 36
               movddup 744(pB), rB1
            #elif KB == 36
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 36
               movaps 1600(pA), rA0
            #elif KB == 36
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 36
               movaps 1616(pA), rA1
            #elif KB == 36
               prefA(-128(pfA))
            #endif

            mulpd 1584(pA), rB2
            addpd rB2, rC22
            #if KB > 36
               movddup 752(pB), rB2
            #elif KB == 36
               prefA((pfA))
            #endif

         #endif
         #if KB > 36
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1632(pA), rB0
            addpd rB0, rC20
            #if KB > 37
               movddup 760(pB), rB0
            #elif KB == 37
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1632(pA), rB1
            addpd rB1, rC21
            #if KB > 37
               movddup 768(pB), rB1
            #elif KB == 37
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 37
               movaps 1648(pA), rA0
            #elif KB == 37
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 37
               movaps 1664(pA), rA1
            #elif KB == 37
               prefA(-128(pfA))
            #endif

            mulpd 1632(pA), rB2
            addpd rB2, rC22
            #if KB > 37
               movddup 776(pB), rB2
            #elif KB == 37
               prefA((pfA))
            #endif

         #endif
         #if KB > 37
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1680(pA), rB0
            addpd rB0, rC20
            #if KB > 38
               movddup 784(pB), rB0
            #elif KB == 38
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1680(pA), rB1
            addpd rB1, rC21
            #if KB > 38
               movddup 792(pB), rB1
            #elif KB == 38
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 38
               movaps 1696(pA), rA0
            #elif KB == 38
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 38
               movaps 1712(pA), rA1
            #elif KB == 38
               prefA(-128(pfA))
            #endif

            mulpd 1680(pA), rB2
            addpd rB2, rC22
            #if KB > 38
               movddup 800(pB), rB2
            #elif KB == 38
               prefA((pfA))
            #endif

         #endif
         #if KB > 38
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1728(pA), rB0
            addpd rB0, rC20
            #if KB > 39
               movddup 808(pB), rB0
            #elif KB == 39
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1728(pA), rB1
            addpd rB1, rC21
            #if KB > 39
               movddup 816(pB), rB1
            #elif KB == 39
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 39
               movaps 1744(pA), rA0
            #elif KB == 39
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 39
               movaps 1760(pA), rA1
            #elif KB == 39
               prefA(-128(pfA))
            #endif

            mulpd 1728(pA), rB2
            addpd rB2, rC22
            #if KB > 39
               movddup 824(pB), rB2
            #elif KB == 39
               prefA((pfA))
            #endif

         #endif
         #if KB > 39
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1776(pA), rB0
            addpd rB0, rC20
            #if KB > 40
               movddup 832(pB), rB0
            #elif KB == 40
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1776(pA), rB1
            addpd rB1, rC21
            #if KB > 40
               movddup 840(pB), rB1
            #elif KB == 40
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 40
               movaps 1792(pA), rA0
            #elif KB == 40
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 40
               movaps 1808(pA), rA1
            #elif KB == 40
               prefA(-128(pfA))
            #endif

            mulpd 1776(pA), rB2
            addpd rB2, rC22
            #if KB > 40
               movddup 848(pB), rB2
            #elif KB == 40
               prefA((pfA))
            #endif

         #endif
         #if KB > 40
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1824(pA), rB0
            addpd rB0, rC20
            #if KB > 41
               movddup 856(pB), rB0
            #elif KB == 41
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1824(pA), rB1
            addpd rB1, rC21
            #if KB > 41
               movddup 864(pB), rB1
            #elif KB == 41
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 41
               movaps 1840(pA), rA0
            #elif KB == 41
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 41
               movaps 1856(pA), rA1
            #elif KB == 41
               prefA(-128(pfA))
            #endif

            mulpd 1824(pA), rB2
            addpd rB2, rC22
            #if KB > 41
               movddup 872(pB), rB2
            #elif KB == 41
               prefA((pfA))
            #endif

         #endif
         #if KB > 41
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1872(pA), rB0
            addpd rB0, rC20
            #if KB > 42
               movddup 880(pB), rB0
            #elif KB == 42
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1872(pA), rB1
            addpd rB1, rC21
            #if KB > 42
               movddup 888(pB), rB1
            #elif KB == 42
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 42
               movaps 1888(pA), rA0
            #elif KB == 42
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 42
               movaps 1904(pA), rA1
            #elif KB == 42
               prefA(-128(pfA))
            #endif

            mulpd 1872(pA), rB2
            addpd rB2, rC22
            #if KB > 42
               movddup 896(pB), rB2
            #elif KB == 42
               prefA((pfA))
            #endif

         #endif
         #if KB > 42
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1920(pA), rB0
            addpd rB0, rC20
            #if KB > 43
               movddup 904(pB), rB0
            #elif KB == 43
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1920(pA), rB1
            addpd rB1, rC21
            #if KB > 43
               movddup 912(pB), rB1
            #elif KB == 43
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 43
               movaps 1936(pA), rA0
            #elif KB == 43
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 43
               movaps 1952(pA), rA1
            #elif KB == 43
               prefA(-128(pfA))
            #endif

            mulpd 1920(pA), rB2
            addpd rB2, rC22
            #if KB > 43
               movddup 920(pB), rB2
            #elif KB == 43
               prefA((pfA))
            #endif

         #endif
         #if KB > 43
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 1968(pA), rB0
            addpd rB0, rC20
            #if KB > 44
               movddup 928(pB), rB0
            #elif KB == 44
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 1968(pA), rB1
            addpd rB1, rC21
            #if KB > 44
               movddup 936(pB), rB1
            #elif KB == 44
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 44
               movaps 1984(pA), rA0
            #elif KB == 44
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 44
               movaps 2000(pA), rA1
            #elif KB == 44
               prefA(-128(pfA))
            #endif

            mulpd 1968(pA), rB2
            addpd rB2, rC22
            #if KB > 44
               movddup 944(pB), rB2
            #elif KB == 44
               prefA((pfA))
            #endif

         #endif
         #if KB > 44
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2016(pA), rB0
            addpd rB0, rC20
            #if KB > 45
               movddup 952(pB), rB0
            #elif KB == 45
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2016(pA), rB1
            addpd rB1, rC21
            #if KB > 45
               movddup 960(pB), rB1
            #elif KB == 45
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 45
               movaps 2032(pA), rA0
            #elif KB == 45
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 45
               movaps 2048(pA), rA1
            #elif KB == 45
               prefA(-128(pfA))
            #endif

            mulpd 2016(pA), rB2
            addpd rB2, rC22
            #if KB > 45
               movddup 968(pB), rB2
            #elif KB == 45
               prefA((pfA))
            #endif

         #endif
         #if KB > 45
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2064(pA), rB0
            addpd rB0, rC20
            #if KB > 46
               movddup 976(pB), rB0
            #elif KB == 46
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2064(pA), rB1
            addpd rB1, rC21
            #if KB > 46
               movddup 984(pB), rB1
            #elif KB == 46
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 46
               movaps 2080(pA), rA0
            #elif KB == 46
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 46
               movaps 2096(pA), rA1
            #elif KB == 46
               prefA(-128(pfA))
            #endif

            mulpd 2064(pA), rB2
            addpd rB2, rC22
            #if KB > 46
               movddup 992(pB), rB2
            #elif KB == 46
               prefA((pfA))
            #endif

         #endif
         #if KB > 46
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2112(pA), rB0
            addpd rB0, rC20
            #if KB > 47
               movddup 1000(pB), rB0
            #elif KB == 47
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2112(pA), rB1
            addpd rB1, rC21
            #if KB > 47
               movddup 1008(pB), rB1
            #elif KB == 47
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 47
               movaps 2128(pA), rA0
            #elif KB == 47
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 47
               movaps 2144(pA), rA1
            #elif KB == 47
               prefA(-128(pfA))
            #endif

            mulpd 2112(pA), rB2
            addpd rB2, rC22
            #if KB > 47
               movddup 1016(pB), rB2
            #elif KB == 47
               prefA((pfA))
            #endif

         #endif
         #if KB > 47
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2160(pA), rB0
            addpd rB0, rC20
            #if KB > 48
               movddup 1024(pB), rB0
            #elif KB == 48
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2160(pA), rB1
            addpd rB1, rC21
            #if KB > 48
               movddup 1032(pB), rB1
            #elif KB == 48
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 48
               movaps 2176(pA), rA0
            #elif KB == 48
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 48
               movaps 2192(pA), rA1
            #elif KB == 48
               prefA(-128(pfA))
            #endif

            mulpd 2160(pA), rB2
            addpd rB2, rC22
            #if KB > 48
               movddup 1040(pB), rB2
            #elif KB == 48
               prefA((pfA))
            #endif

         #endif
         #if KB > 48
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2208(pA), rB0
            addpd rB0, rC20
            #if KB > 49
               movddup 1048(pB), rB0
            #elif KB == 49
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2208(pA), rB1
            addpd rB1, rC21
            #if KB > 49
               movddup 1056(pB), rB1
            #elif KB == 49
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 49
               movaps 2224(pA), rA0
            #elif KB == 49
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 49
               movaps 2240(pA), rA1
            #elif KB == 49
               prefA(-128(pfA))
            #endif

            mulpd 2208(pA), rB2
            addpd rB2, rC22
            #if KB > 49
               movddup 1064(pB), rB2
            #elif KB == 49
               prefA((pfA))
            #endif

         #endif
         #if KB > 49
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2256(pA), rB0
            addpd rB0, rC20
            #if KB > 50
               movddup 1072(pB), rB0
            #elif KB == 50
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2256(pA), rB1
            addpd rB1, rC21
            #if KB > 50
               movddup 1080(pB), rB1
            #elif KB == 50
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 50
               movaps 2272(pA), rA0
            #elif KB == 50
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 50
               movaps 2288(pA), rA1
            #elif KB == 50
               prefA(-128(pfA))
            #endif

            mulpd 2256(pA), rB2
            addpd rB2, rC22
            #if KB > 50
               movddup 1088(pB), rB2
            #elif KB == 50
               prefA((pfA))
            #endif

         #endif
         #if KB > 50
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2304(pA), rB0
            addpd rB0, rC20
            #if KB > 51
               movddup 1096(pB), rB0
            #elif KB == 51
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2304(pA), rB1
            addpd rB1, rC21
            #if KB > 51
               movddup 1104(pB), rB1
            #elif KB == 51
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 51
               movaps 2320(pA), rA0
            #elif KB == 51
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 51
               movaps 2336(pA), rA1
            #elif KB == 51
               prefA(-128(pfA))
            #endif

            mulpd 2304(pA), rB2
            addpd rB2, rC22
            #if KB > 51
               movddup 1112(pB), rB2
            #elif KB == 51
               prefA((pfA))
            #endif

         #endif
         #if KB > 51
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2352(pA), rB0
            addpd rB0, rC20
            #if KB > 52
               movddup 1120(pB), rB0
            #elif KB == 52
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2352(pA), rB1
            addpd rB1, rC21
            #if KB > 52
               movddup 1128(pB), rB1
            #elif KB == 52
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 52
               movaps 2368(pA), rA0
            #elif KB == 52
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 52
               movaps 2384(pA), rA1
            #elif KB == 52
               prefA(-128(pfA))
            #endif

            mulpd 2352(pA), rB2
            addpd rB2, rC22
            #if KB > 52
               movddup 1136(pB), rB2
            #elif KB == 52
               prefA((pfA))
            #endif

         #endif
         #if KB > 52
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2400(pA), rB0
            addpd rB0, rC20
            #if KB > 53
               movddup 1144(pB), rB0
            #elif KB == 53
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2400(pA), rB1
            addpd rB1, rC21
            #if KB > 53
               movddup 1152(pB), rB1
            #elif KB == 53
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 53
               movaps 2416(pA), rA0
            #elif KB == 53
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 53
               movaps 2432(pA), rA1
            #elif KB == 53
               prefA(-128(pfA))
            #endif

            mulpd 2400(pA), rB2
            addpd rB2, rC22
            #if KB > 53
               movddup 1160(pB), rB2
            #elif KB == 53
               prefA((pfA))
            #endif

         #endif
         #if KB > 53
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2448(pA), rB0
            addpd rB0, rC20
            #if KB > 54
               movddup 1168(pB), rB0
            #elif KB == 54
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2448(pA), rB1
            addpd rB1, rC21
            #if KB > 54
               movddup 1176(pB), rB1
            #elif KB == 54
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 54
               movaps 2464(pA), rA0
            #elif KB == 54
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 54
               movaps 2480(pA), rA1
            #elif KB == 54
               prefA(-128(pfA))
            #endif

            mulpd 2448(pA), rB2
            addpd rB2, rC22
            #if KB > 54
               movddup 1184(pB), rB2
            #elif KB == 54
               prefA((pfA))
            #endif

         #endif
         #if KB > 54
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2496(pA), rB0
            addpd rB0, rC20
            #if KB > 55
               movddup 1192(pB), rB0
            #elif KB == 55
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2496(pA), rB1
            addpd rB1, rC21
            #if KB > 55
               movddup 1200(pB), rB1
            #elif KB == 55
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 55
               movaps 2512(pA), rA0
            #elif KB == 55
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 55
               movaps 2528(pA), rA1
            #elif KB == 55
               prefA(-128(pfA))
            #endif

            mulpd 2496(pA), rB2
            addpd rB2, rC22
            #if KB > 55
               movddup 1208(pB), rB2
            #elif KB == 55
               prefA((pfA))
            #endif

         #endif
         #if KB > 55
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2544(pA), rB0
            addpd rB0, rC20
            #if KB > 56
               movddup 1216(pB), rB0
            #elif KB == 56
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2544(pA), rB1
            addpd rB1, rC21
            #if KB > 56
               movddup 1224(pB), rB1
            #elif KB == 56
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 56
               movaps 2560(pA), rA0
            #elif KB == 56
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 56
               movaps 2576(pA), rA1
            #elif KB == 56
               prefA(-128(pfA))
            #endif

            mulpd 2544(pA), rB2
            addpd rB2, rC22
            #if KB > 56
               movddup 1232(pB), rB2
            #elif KB == 56
               prefA((pfA))
            #endif

         #endif
         #if KB > 56
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2592(pA), rB0
            addpd rB0, rC20
            #if KB > 57
               movddup 1240(pB), rB0
            #elif KB == 57
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2592(pA), rB1
            addpd rB1, rC21
            #if KB > 57
               movddup 1248(pB), rB1
            #elif KB == 57
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 57
               movaps 2608(pA), rA0
            #elif KB == 57
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 57
               movaps 2624(pA), rA1
            #elif KB == 57
               prefA(-128(pfA))
            #endif

            mulpd 2592(pA), rB2
            addpd rB2, rC22
            #if KB > 57
               movddup 1256(pB), rB2
            #elif KB == 57
               prefA((pfA))
            #endif

         #endif
         #if KB > 57
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2640(pA), rB0
            addpd rB0, rC20
            #if KB > 58
               movddup 1264(pB), rB0
            #elif KB == 58
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2640(pA), rB1
            addpd rB1, rC21
            #if KB > 58
               movddup 1272(pB), rB1
            #elif KB == 58
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 58
               movaps 2656(pA), rA0
            #elif KB == 58
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 58
               movaps 2672(pA), rA1
            #elif KB == 58
               prefA(-128(pfA))
            #endif

            mulpd 2640(pA), rB2
            addpd rB2, rC22
            #if KB > 58
               movddup 1280(pB), rB2
            #elif KB == 58
               prefA((pfA))
            #endif

         #endif
         #if KB > 58
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2688(pA), rB0
            addpd rB0, rC20
            #if KB > 59
               movddup 1288(pB), rB0
            #elif KB == 59
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2688(pA), rB1
            addpd rB1, rC21
            #if KB > 59
               movddup 1296(pB), rB1
            #elif KB == 59
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 59
               movaps 2704(pA), rA0
            #elif KB == 59
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 59
               movaps 2720(pA), rA1
            #elif KB == 59
               prefA(-128(pfA))
            #endif

            mulpd 2688(pA), rB2
            addpd rB2, rC22
            #if KB > 59
               movddup 1304(pB), rB2
            #elif KB == 59
               prefA((pfA))
            #endif

         #endif
         #if KB > 59
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2736(pA), rB0
            addpd rB0, rC20
            #if KB > 60
               movddup 1312(pB), rB0
            #elif KB == 60
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2736(pA), rB1
            addpd rB1, rC21
            #if KB > 60
               movddup 1320(pB), rB1
            #elif KB == 60
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 60
               movaps 2752(pA), rA0
            #elif KB == 60
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 60
               movaps 2768(pA), rA1
            #elif KB == 60
               prefA(-128(pfA))
            #endif

            mulpd 2736(pA), rB2
            addpd rB2, rC22
            #if KB > 60
               movddup 1328(pB), rB2
            #elif KB == 60
               prefA((pfA))
            #endif

         #endif
         #if KB > 60
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2784(pA), rB0
            addpd rB0, rC20
            #if KB > 61
               movddup 1336(pB), rB0
            #elif KB == 61
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2784(pA), rB1
            addpd rB1, rC21
            #if KB > 61
               movddup 1344(pB), rB1
            #elif KB == 61
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 61
               movaps 2800(pA), rA0
            #elif KB == 61
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 61
               movaps 2816(pA), rA1
            #elif KB == 61
               prefA(-128(pfA))
            #endif

            mulpd 2784(pA), rB2
            addpd rB2, rC22
            #if KB > 61
               movddup 1352(pB), rB2
            #elif KB == 61
               prefA((pfA))
            #endif

         #endif
         #if KB > 61
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2832(pA), rB0
            addpd rB0, rC20
            #if KB > 62
               movddup 1360(pB), rB0
            #elif KB == 62
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2832(pA), rB1
            addpd rB1, rC21
            #if KB > 62
               movddup 1368(pB), rB1
            #elif KB == 62
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 62
               movaps 2848(pA), rA0
            #elif KB == 62
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 62
               movaps 2864(pA), rA1
            #elif KB == 62
               prefA(-128(pfA))
            #endif

            mulpd 2832(pA), rB2
            addpd rB2, rC22
            #if KB > 62
               movddup 1376(pB), rB2
            #elif KB == 62
               prefA((pfA))
            #endif

         #endif
         #if KB > 62
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2880(pA), rB0
            addpd rB0, rC20
            #if KB > 63
               movddup 1384(pB), rB0
            #elif KB == 63
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2880(pA), rB1
            addpd rB1, rC21
            #if KB > 63
               movddup 1392(pB), rB1
            #elif KB == 63
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 63
               movaps 2896(pA), rA0
            #elif KB == 63
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 63
               movaps 2912(pA), rA1
            #elif KB == 63
               prefA(-128(pfA))
            #endif

            mulpd 2880(pA), rB2
            addpd rB2, rC22
            #if KB > 63
               movddup 1400(pB), rB2
            #elif KB == 63
               prefA((pfA))
            #endif

         #endif
         #if KB > 63
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2928(pA), rB0
            addpd rB0, rC20
            #if KB > 64
               movddup 1408(pB), rB0
            #elif KB == 64
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2928(pA), rB1
            addpd rB1, rC21
            #if KB > 64
               movddup 1416(pB), rB1
            #elif KB == 64
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 64
               movaps 2944(pA), rA0
            #elif KB == 64
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 64
               movaps 2960(pA), rA1
            #elif KB == 64
               prefA(-128(pfA))
            #endif

            mulpd 2928(pA), rB2
            addpd rB2, rC22
            #if KB > 64
               movddup 1424(pB), rB2
            #elif KB == 64
               prefA((pfA))
            #endif

         #endif
         #if KB > 64
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 2976(pA), rB0
            addpd rB0, rC20
            #if KB > 65
               movddup 1432(pB), rB0
            #elif KB == 65
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 2976(pA), rB1
            addpd rB1, rC21
            #if KB > 65
               movddup 1440(pB), rB1
            #elif KB == 65
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 65
               movaps 2992(pA), rA0
            #elif KB == 65
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 65
               movaps 3008(pA), rA1
            #elif KB == 65
               prefA(-128(pfA))
            #endif

            mulpd 2976(pA), rB2
            addpd rB2, rC22
            #if KB > 65
               movddup 1448(pB), rB2
            #elif KB == 65
               prefA((pfA))
            #endif

         #endif
         #if KB > 65
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3024(pA), rB0
            addpd rB0, rC20
            #if KB > 66
               movddup 1456(pB), rB0
            #elif KB == 66
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3024(pA), rB1
            addpd rB1, rC21
            #if KB > 66
               movddup 1464(pB), rB1
            #elif KB == 66
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 66
               movaps 3040(pA), rA0
            #elif KB == 66
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 66
               movaps 3056(pA), rA1
            #elif KB == 66
               prefA(-128(pfA))
            #endif

            mulpd 3024(pA), rB2
            addpd rB2, rC22
            #if KB > 66
               movddup 1472(pB), rB2
            #elif KB == 66
               prefA((pfA))
            #endif

         #endif
         #if KB > 66
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3072(pA), rB0
            addpd rB0, rC20
            #if KB > 67
               movddup 1480(pB), rB0
            #elif KB == 67
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3072(pA), rB1
            addpd rB1, rC21
            #if KB > 67
               movddup 1488(pB), rB1
            #elif KB == 67
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 67
               movaps 3088(pA), rA0
            #elif KB == 67
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 67
               movaps 3104(pA), rA1
            #elif KB == 67
               prefA(-128(pfA))
            #endif

            mulpd 3072(pA), rB2
            addpd rB2, rC22
            #if KB > 67
               movddup 1496(pB), rB2
            #elif KB == 67
               prefA((pfA))
            #endif

         #endif
         #if KB > 67
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3120(pA), rB0
            addpd rB0, rC20
            #if KB > 68
               movddup 1504(pB), rB0
            #elif KB == 68
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3120(pA), rB1
            addpd rB1, rC21
            #if KB > 68
               movddup 1512(pB), rB1
            #elif KB == 68
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 68
               movaps 3136(pA), rA0
            #elif KB == 68
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 68
               movaps 3152(pA), rA1
            #elif KB == 68
               prefA(-128(pfA))
            #endif

            mulpd 3120(pA), rB2
            addpd rB2, rC22
            #if KB > 68
               movddup 1520(pB), rB2
            #elif KB == 68
               prefA((pfA))
            #endif

         #endif
         #if KB > 68
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3168(pA), rB0
            addpd rB0, rC20
            #if KB > 69
               movddup 1528(pB), rB0
            #elif KB == 69
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3168(pA), rB1
            addpd rB1, rC21
            #if KB > 69
               movddup 1536(pB), rB1
            #elif KB == 69
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 69
               movaps 3184(pA), rA0
            #elif KB == 69
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 69
               movaps 3200(pA), rA1
            #elif KB == 69
               prefA(-128(pfA))
            #endif

            mulpd 3168(pA), rB2
            addpd rB2, rC22
            #if KB > 69
               movddup 1544(pB), rB2
            #elif KB == 69
               prefA((pfA))
            #endif

         #endif
         #if KB > 69
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3216(pA), rB0
            addpd rB0, rC20
            #if KB > 70
               movddup 1552(pB), rB0
            #elif KB == 70
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3216(pA), rB1
            addpd rB1, rC21
            #if KB > 70
               movddup 1560(pB), rB1
            #elif KB == 70
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 70
               movaps 3232(pA), rA0
            #elif KB == 70
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 70
               movaps 3248(pA), rA1
            #elif KB == 70
               prefA(-128(pfA))
            #endif

            mulpd 3216(pA), rB2
            addpd rB2, rC22
            #if KB > 70
               movddup 1568(pB), rB2
            #elif KB == 70
               prefA((pfA))
            #endif

         #endif
         #if KB > 70
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3264(pA), rB0
            addpd rB0, rC20
            #if KB > 71
               movddup 1576(pB), rB0
            #elif KB == 71
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3264(pA), rB1
            addpd rB1, rC21
            #if KB > 71
               movddup 1584(pB), rB1
            #elif KB == 71
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 71
               movaps 3280(pA), rA0
            #elif KB == 71
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 71
               movaps 3296(pA), rA1
            #elif KB == 71
               prefA(-128(pfA))
            #endif

            mulpd 3264(pA), rB2
            addpd rB2, rC22
            #if KB > 71
               movddup 1592(pB), rB2
            #elif KB == 71
               prefA((pfA))
            #endif

         #endif
         #if KB > 71
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3312(pA), rB0
            addpd rB0, rC20
            #if KB > 72
               movddup 1600(pB), rB0
            #elif KB == 72
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3312(pA), rB1
            addpd rB1, rC21
            #if KB > 72
               movddup 1608(pB), rB1
            #elif KB == 72
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 72
               movaps 3328(pA), rA0
            #elif KB == 72
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 72
               movaps 3344(pA), rA1
            #elif KB == 72
               prefA(-128(pfA))
            #endif

            mulpd 3312(pA), rB2
            addpd rB2, rC22
            #if KB > 72
               movddup 1616(pB), rB2
            #elif KB == 72
               prefA((pfA))
            #endif

         #endif
         #if KB > 72
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3360(pA), rB0
            addpd rB0, rC20
            #if KB > 73
               movddup 1624(pB), rB0
            #elif KB == 73
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3360(pA), rB1
            addpd rB1, rC21
            #if KB > 73
               movddup 1632(pB), rB1
            #elif KB == 73
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 73
               movaps 3376(pA), rA0
            #elif KB == 73
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 73
               movaps 3392(pA), rA1
            #elif KB == 73
               prefA(-128(pfA))
            #endif

            mulpd 3360(pA), rB2
            addpd rB2, rC22
            #if KB > 73
               movddup 1640(pB), rB2
            #elif KB == 73
               prefA((pfA))
            #endif

         #endif
         #if KB > 73
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3408(pA), rB0
            addpd rB0, rC20
            #if KB > 74
               movddup 1648(pB), rB0
            #elif KB == 74
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3408(pA), rB1
            addpd rB1, rC21
            #if KB > 74
               movddup 1656(pB), rB1
            #elif KB == 74
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 74
               movaps 3424(pA), rA0
            #elif KB == 74
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 74
               movaps 3440(pA), rA1
            #elif KB == 74
               prefA(-128(pfA))
            #endif

            mulpd 3408(pA), rB2
            addpd rB2, rC22
            #if KB > 74
               movddup 1664(pB), rB2
            #elif KB == 74
               prefA((pfA))
            #endif

         #endif
         #if KB > 74
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3456(pA), rB0
            addpd rB0, rC20
            #if KB > 75
               movddup 1672(pB), rB0
            #elif KB == 75
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3456(pA), rB1
            addpd rB1, rC21
            #if KB > 75
               movddup 1680(pB), rB1
            #elif KB == 75
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 75
               movaps 3472(pA), rA0
            #elif KB == 75
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 75
               movaps 3488(pA), rA1
            #elif KB == 75
               prefA(-128(pfA))
            #endif

            mulpd 3456(pA), rB2
            addpd rB2, rC22
            #if KB > 75
               movddup 1688(pB), rB2
            #elif KB == 75
               prefA((pfA))
            #endif

         #endif
         #if KB > 75
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3504(pA), rB0
            addpd rB0, rC20
            #if KB > 76
               movddup 1696(pB), rB0
            #elif KB == 76
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3504(pA), rB1
            addpd rB1, rC21
            #if KB > 76
               movddup 1704(pB), rB1
            #elif KB == 76
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 76
               movaps 3520(pA), rA0
            #elif KB == 76
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 76
               movaps 3536(pA), rA1
            #elif KB == 76
               prefA(-128(pfA))
            #endif

            mulpd 3504(pA), rB2
            addpd rB2, rC22
            #if KB > 76
               movddup 1712(pB), rB2
            #elif KB == 76
               prefA((pfA))
            #endif

         #endif
         #if KB > 76
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3552(pA), rB0
            addpd rB0, rC20
            #if KB > 77
               movddup 1720(pB), rB0
            #elif KB == 77
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3552(pA), rB1
            addpd rB1, rC21
            #if KB > 77
               movddup 1728(pB), rB1
            #elif KB == 77
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 77
               movaps 3568(pA), rA0
            #elif KB == 77
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 77
               movaps 3584(pA), rA1
            #elif KB == 77
               prefA(-128(pfA))
            #endif

            mulpd 3552(pA), rB2
            addpd rB2, rC22
            #if KB > 77
               movddup 1736(pB), rB2
            #elif KB == 77
               prefA((pfA))
            #endif

         #endif
         #if KB > 77
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3600(pA), rB0
            addpd rB0, rC20
            #if KB > 78
               movddup 1744(pB), rB0
            #elif KB == 78
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3600(pA), rB1
            addpd rB1, rC21
            #if KB > 78
               movddup 1752(pB), rB1
            #elif KB == 78
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 78
               movaps 3616(pA), rA0
            #elif KB == 78
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 78
               movaps 3632(pA), rA1
            #elif KB == 78
               prefA(-128(pfA))
            #endif

            mulpd 3600(pA), rB2
            addpd rB2, rC22
            #if KB > 78
               movddup 1760(pB), rB2
            #elif KB == 78
               prefA((pfA))
            #endif

         #endif
         #if KB > 78
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3648(pA), rB0
            addpd rB0, rC20
            #if KB > 79
               movddup 1768(pB), rB0
            #elif KB == 79
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3648(pA), rB1
            addpd rB1, rC21
            #if KB > 79
               movddup 1776(pB), rB1
            #elif KB == 79
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 79
               movaps 3664(pA), rA0
            #elif KB == 79
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 79
               movaps 3680(pA), rA1
            #elif KB == 79
               prefA(-128(pfA))
            #endif

            mulpd 3648(pA), rB2
            addpd rB2, rC22
            #if KB > 79
               movddup 1784(pB), rB2
            #elif KB == 79
               prefA((pfA))
            #endif

         #endif
         #if KB > 79
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3696(pA), rB0
            addpd rB0, rC20
            #if KB > 80
               movddup 1792(pB), rB0
            #elif KB == 80
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3696(pA), rB1
            addpd rB1, rC21
            #if KB > 80
               movddup 1800(pB), rB1
            #elif KB == 80
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 80
               movaps 3712(pA), rA0
            #elif KB == 80
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 80
               movaps 3728(pA), rA1
            #elif KB == 80
               prefA(-128(pfA))
            #endif

            mulpd 3696(pA), rB2
            addpd rB2, rC22
            #if KB > 80
               movddup 1808(pB), rB2
            #elif KB == 80
               prefA((pfA))
            #endif

         #endif
         #if KB > 80
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3744(pA), rB0
            addpd rB0, rC20
            #if KB > 81
               movddup 1816(pB), rB0
            #elif KB == 81
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3744(pA), rB1
            addpd rB1, rC21
            #if KB > 81
               movddup 1824(pB), rB1
            #elif KB == 81
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 81
               movaps 3760(pA), rA0
            #elif KB == 81
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 81
               movaps 3776(pA), rA1
            #elif KB == 81
               prefA(-128(pfA))
            #endif

            mulpd 3744(pA), rB2
            addpd rB2, rC22
            #if KB > 81
               movddup 1832(pB), rB2
            #elif KB == 81
               prefA((pfA))
            #endif

         #endif
         #if KB > 81
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3792(pA), rB0
            addpd rB0, rC20
            #if KB > 82
               movddup 1840(pB), rB0
            #elif KB == 82
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3792(pA), rB1
            addpd rB1, rC21
            #if KB > 82
               movddup 1848(pB), rB1
            #elif KB == 82
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 82
               movaps 3808(pA), rA0
            #elif KB == 82
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 82
               movaps 3824(pA), rA1
            #elif KB == 82
               prefA(-128(pfA))
            #endif

            mulpd 3792(pA), rB2
            addpd rB2, rC22
            #if KB > 82
               movddup 1856(pB), rB2
            #elif KB == 82
               prefA((pfA))
            #endif

         #endif
         #if KB > 82
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3840(pA), rB0
            addpd rB0, rC20
            #if KB > 83
               movddup 1864(pB), rB0
            #elif KB == 83
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3840(pA), rB1
            addpd rB1, rC21
            #if KB > 83
               movddup 1872(pB), rB1
            #elif KB == 83
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 83
               movaps 3856(pA), rA0
            #elif KB == 83
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 83
               movaps 3872(pA), rA1
            #elif KB == 83
               prefA(-128(pfA))
            #endif

            mulpd 3840(pA), rB2
            addpd rB2, rC22
            #if KB > 83
               movddup 1880(pB), rB2
            #elif KB == 83
               prefA((pfA))
            #endif

         #endif
         #if KB > 83
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3888(pA), rB0
            addpd rB0, rC20
            #if KB > 84
               movddup 1888(pB), rB0
            #elif KB == 84
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3888(pA), rB1
            addpd rB1, rC21
            #if KB > 84
               movddup 1896(pB), rB1
            #elif KB == 84
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 84
               movaps 3904(pA), rA0
            #elif KB == 84
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 84
               movaps 3920(pA), rA1
            #elif KB == 84
               prefA(-128(pfA))
            #endif

            mulpd 3888(pA), rB2
            addpd rB2, rC22
            #if KB > 84
               movddup 1904(pB), rB2
            #elif KB == 84
               prefA((pfA))
            #endif

         #endif
         #if KB > 84
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3936(pA), rB0
            addpd rB0, rC20
            #if KB > 85
               movddup 1912(pB), rB0
            #elif KB == 85
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3936(pA), rB1
            addpd rB1, rC21
            #if KB > 85
               movddup 1920(pB), rB1
            #elif KB == 85
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 85
               movaps 3952(pA), rA0
            #elif KB == 85
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 85
               movaps 3968(pA), rA1
            #elif KB == 85
               prefA(-128(pfA))
            #endif

            mulpd 3936(pA), rB2
            addpd rB2, rC22
            #if KB > 85
               movddup 1928(pB), rB2
            #elif KB == 85
               prefA((pfA))
            #endif

         #endif
         #if KB > 85
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 3984(pA), rB0
            addpd rB0, rC20
            #if KB > 86
               movddup 1936(pB), rB0
            #elif KB == 86
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 3984(pA), rB1
            addpd rB1, rC21
            #if KB > 86
               movddup 1944(pB), rB1
            #elif KB == 86
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 86
               movaps 4000(pA), rA0
            #elif KB == 86
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 86
               movaps 4016(pA), rA1
            #elif KB == 86
               prefA(-128(pfA))
            #endif

            mulpd 3984(pA), rB2
            addpd rB2, rC22
            #if KB > 86
               movddup 1952(pB), rB2
            #elif KB == 86
               prefA((pfA))
            #endif

         #endif
         #if KB > 86
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4032(pA), rB0
            addpd rB0, rC20
            #if KB > 87
               movddup 1960(pB), rB0
            #elif KB == 87
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4032(pA), rB1
            addpd rB1, rC21
            #if KB > 87
               movddup 1968(pB), rB1
            #elif KB == 87
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 87
               movaps 4048(pA), rA0
            #elif KB == 87
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 87
               movaps 4064(pA), rA1
            #elif KB == 87
               prefA(-128(pfA))
            #endif

            mulpd 4032(pA), rB2
            addpd rB2, rC22
            #if KB > 87
               movddup 1976(pB), rB2
            #elif KB == 87
               prefA((pfA))
            #endif

         #endif
         #if KB > 87
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4080(pA), rB0
            addpd rB0, rC20
            #if KB > 88
               movddup 1984(pB), rB0
            #elif KB == 88
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4080(pA), rB1
            addpd rB1, rC21
            #if KB > 88
               movddup 1992(pB), rB1
            #elif KB == 88
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 88
               movaps 4096(pA), rA0
            #elif KB == 88
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 88
               movaps 4112(pA), rA1
            #elif KB == 88
               prefA(-128(pfA))
            #endif

            mulpd 4080(pA), rB2
            addpd rB2, rC22
            #if KB > 88
               movddup 2000(pB), rB2
            #elif KB == 88
               prefA((pfA))
            #endif

         #endif
         #if KB > 88
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4128(pA), rB0
            addpd rB0, rC20
            #if KB > 89
               movddup 2008(pB), rB0
            #elif KB == 89
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4128(pA), rB1
            addpd rB1, rC21
            #if KB > 89
               movddup 2016(pB), rB1
            #elif KB == 89
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 89
               movaps 4144(pA), rA0
            #elif KB == 89
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 89
               movaps 4160(pA), rA1
            #elif KB == 89
               prefA(-128(pfA))
            #endif

            mulpd 4128(pA), rB2
            addpd rB2, rC22
            #if KB > 89
               movddup 2024(pB), rB2
            #elif KB == 89
               prefA((pfA))
            #endif

         #endif
         #if KB > 89
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4176(pA), rB0
            addpd rB0, rC20
            #if KB > 90
               movddup 2032(pB), rB0
            #elif KB == 90
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4176(pA), rB1
            addpd rB1, rC21
            #if KB > 90
               movddup 2040(pB), rB1
            #elif KB == 90
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 90
               movaps 4192(pA), rA0
            #elif KB == 90
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 90
               movaps 4208(pA), rA1
            #elif KB == 90
               prefA(-128(pfA))
            #endif

            mulpd 4176(pA), rB2
            addpd rB2, rC22
            #if KB > 90
               movddup 2048(pB), rB2
            #elif KB == 90
               prefA((pfA))
            #endif

         #endif
         #if KB > 90
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4224(pA), rB0
            addpd rB0, rC20
            #if KB > 91
               movddup 2056(pB), rB0
            #elif KB == 91
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4224(pA), rB1
            addpd rB1, rC21
            #if KB > 91
               movddup 2064(pB), rB1
            #elif KB == 91
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 91
               movaps 4240(pA), rA0
            #elif KB == 91
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 91
               movaps 4256(pA), rA1
            #elif KB == 91
               prefA(-128(pfA))
            #endif

            mulpd 4224(pA), rB2
            addpd rB2, rC22
            #if KB > 91
               movddup 2072(pB), rB2
            #elif KB == 91
               prefA((pfA))
            #endif

         #endif
         #if KB > 91
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4272(pA), rB0
            addpd rB0, rC20
            #if KB > 92
               movddup 2080(pB), rB0
            #elif KB == 92
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4272(pA), rB1
            addpd rB1, rC21
            #if KB > 92
               movddup 2088(pB), rB1
            #elif KB == 92
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 92
               movaps 4288(pA), rA0
            #elif KB == 92
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 92
               movaps 4304(pA), rA1
            #elif KB == 92
               prefA(-128(pfA))
            #endif

            mulpd 4272(pA), rB2
            addpd rB2, rC22
            #if KB > 92
               movddup 2096(pB), rB2
            #elif KB == 92
               prefA((pfA))
            #endif

         #endif
         #if KB > 92
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4320(pA), rB0
            addpd rB0, rC20
            #if KB > 93
               movddup 2104(pB), rB0
            #elif KB == 93
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4320(pA), rB1
            addpd rB1, rC21
            #if KB > 93
               movddup 2112(pB), rB1
            #elif KB == 93
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 93
               movaps 4336(pA), rA0
            #elif KB == 93
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 93
               movaps 4352(pA), rA1
            #elif KB == 93
               prefA(-128(pfA))
            #endif

            mulpd 4320(pA), rB2
            addpd rB2, rC22
            #if KB > 93
               movddup 2120(pB), rB2
            #elif KB == 93
               prefA((pfA))
            #endif

         #endif
         #if KB > 93
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4368(pA), rB0
            addpd rB0, rC20
            #if KB > 94
               movddup 2128(pB), rB0
            #elif KB == 94
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4368(pA), rB1
            addpd rB1, rC21
            #if KB > 94
               movddup 2136(pB), rB1
            #elif KB == 94
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 94
               movaps 4384(pA), rA0
            #elif KB == 94
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 94
               movaps 4400(pA), rA1
            #elif KB == 94
               prefA(-128(pfA))
            #endif

            mulpd 4368(pA), rB2
            addpd rB2, rC22
            #if KB > 94
               movddup 2144(pB), rB2
            #elif KB == 94
               prefA((pfA))
            #endif

         #endif
         #if KB > 94
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4416(pA), rB0
            addpd rB0, rC20
            #if KB > 95
               movddup 2152(pB), rB0
            #elif KB == 95
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4416(pA), rB1
            addpd rB1, rC21
            #if KB > 95
               movddup 2160(pB), rB1
            #elif KB == 95
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 95
               movaps 4432(pA), rA0
            #elif KB == 95
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 95
               movaps 4448(pA), rA1
            #elif KB == 95
               prefA(-128(pfA))
            #endif

            mulpd 4416(pA), rB2
            addpd rB2, rC22
            #if KB > 95
               movddup 2168(pB), rB2
            #elif KB == 95
               prefA((pfA))
            #endif

         #endif
         #if KB > 95
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4464(pA), rB0
            addpd rB0, rC20
            #if KB > 96
               movddup 2176(pB), rB0
            #elif KB == 96
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4464(pA), rB1
            addpd rB1, rC21
            #if KB > 96
               movddup 2184(pB), rB1
            #elif KB == 96
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 96
               movaps 4480(pA), rA0
            #elif KB == 96
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 96
               movaps 4496(pA), rA1
            #elif KB == 96
               prefA(-128(pfA))
            #endif

            mulpd 4464(pA), rB2
            addpd rB2, rC22
            #if KB > 96
               movddup 2192(pB), rB2
            #elif KB == 96
               prefA((pfA))
            #endif

         #endif
         #if KB > 96
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4512(pA), rB0
            addpd rB0, rC20
            #if KB > 97
               movddup 2200(pB), rB0
            #elif KB == 97
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4512(pA), rB1
            addpd rB1, rC21
            #if KB > 97
               movddup 2208(pB), rB1
            #elif KB == 97
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 97
               movaps 4528(pA), rA0
            #elif KB == 97
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 97
               movaps 4544(pA), rA1
            #elif KB == 97
               prefA(-128(pfA))
            #endif

            mulpd 4512(pA), rB2
            addpd rB2, rC22
            #if KB > 97
               movddup 2216(pB), rB2
            #elif KB == 97
               prefA((pfA))
            #endif

         #endif
         #if KB > 97
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4560(pA), rB0
            addpd rB0, rC20
            #if KB > 98
               movddup 2224(pB), rB0
            #elif KB == 98
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4560(pA), rB1
            addpd rB1, rC21
            #if KB > 98
               movddup 2232(pB), rB1
            #elif KB == 98
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 98
               movaps 4576(pA), rA0
            #elif KB == 98
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 98
               movaps 4592(pA), rA1
            #elif KB == 98
               prefA(-128(pfA))
            #endif

            mulpd 4560(pA), rB2
            addpd rB2, rC22
            #if KB > 98
               movddup 2240(pB), rB2
            #elif KB == 98
               prefA((pfA))
            #endif

         #endif
         #if KB > 98
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4608(pA), rB0
            addpd rB0, rC20
            #if KB > 99
               movddup 2248(pB), rB0
            #elif KB == 99
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4608(pA), rB1
            addpd rB1, rC21
            #if KB > 99
               movddup 2256(pB), rB1
            #elif KB == 99
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 99
               movaps 4624(pA), rA0
            #elif KB == 99
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 99
               movaps 4640(pA), rA1
            #elif KB == 99
               prefA(-128(pfA))
            #endif

            mulpd 4608(pA), rB2
            addpd rB2, rC22
            #if KB > 99
               movddup 2264(pB), rB2
            #elif KB == 99
               prefA((pfA))
            #endif

         #endif
         #if KB > 99
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4656(pA), rB0
            addpd rB0, rC20
            #if KB > 100
               movddup 2272(pB), rB0
            #elif KB == 100
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4656(pA), rB1
            addpd rB1, rC21
            #if KB > 100
               movddup 2280(pB), rB1
            #elif KB == 100
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 100
               movaps 4672(pA), rA0
            #elif KB == 100
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 100
               movaps 4688(pA), rA1
            #elif KB == 100
               prefA(-128(pfA))
            #endif

            mulpd 4656(pA), rB2
            addpd rB2, rC22
            #if KB > 100
               movddup 2288(pB), rB2
            #elif KB == 100
               prefA((pfA))
            #endif

         #endif
         #if KB > 100
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4704(pA), rB0
            addpd rB0, rC20
            #if KB > 101
               movddup 2296(pB), rB0
            #elif KB == 101
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4704(pA), rB1
            addpd rB1, rC21
            #if KB > 101
               movddup 2304(pB), rB1
            #elif KB == 101
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 101
               movaps 4720(pA), rA0
            #elif KB == 101
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 101
               movaps 4736(pA), rA1
            #elif KB == 101
               prefA(-128(pfA))
            #endif

            mulpd 4704(pA), rB2
            addpd rB2, rC22
            #if KB > 101
               movddup 2312(pB), rB2
            #elif KB == 101
               prefA((pfA))
            #endif

         #endif
         #if KB > 101
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4752(pA), rB0
            addpd rB0, rC20
            #if KB > 102
               movddup 2320(pB), rB0
            #elif KB == 102
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4752(pA), rB1
            addpd rB1, rC21
            #if KB > 102
               movddup 2328(pB), rB1
            #elif KB == 102
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 102
               movaps 4768(pA), rA0
            #elif KB == 102
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 102
               movaps 4784(pA), rA1
            #elif KB == 102
               prefA(-128(pfA))
            #endif

            mulpd 4752(pA), rB2
            addpd rB2, rC22
            #if KB > 102
               movddup 2336(pB), rB2
            #elif KB == 102
               prefA((pfA))
            #endif

         #endif
         #if KB > 102
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4800(pA), rB0
            addpd rB0, rC20
            #if KB > 103
               movddup 2344(pB), rB0
            #elif KB == 103
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4800(pA), rB1
            addpd rB1, rC21
            #if KB > 103
               movddup 2352(pB), rB1
            #elif KB == 103
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 103
               movaps 4816(pA), rA0
            #elif KB == 103
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 103
               movaps 4832(pA), rA1
            #elif KB == 103
               prefA(-128(pfA))
            #endif

            mulpd 4800(pA), rB2
            addpd rB2, rC22
            #if KB > 103
               movddup 2360(pB), rB2
            #elif KB == 103
               prefA((pfA))
            #endif

         #endif
         #if KB > 103
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4848(pA), rB0
            addpd rB0, rC20
            #if KB > 104
               movddup 2368(pB), rB0
            #elif KB == 104
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4848(pA), rB1
            addpd rB1, rC21
            #if KB > 104
               movddup 2376(pB), rB1
            #elif KB == 104
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 104
               movaps 4864(pA), rA0
            #elif KB == 104
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 104
               movaps 4880(pA), rA1
            #elif KB == 104
               prefA(-128(pfA))
            #endif

            mulpd 4848(pA), rB2
            addpd rB2, rC22
            #if KB > 104
               movddup 2384(pB), rB2
            #elif KB == 104
               prefA((pfA))
            #endif

         #endif
         #if KB > 104
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4896(pA), rB0
            addpd rB0, rC20
            #if KB > 105
               movddup 2392(pB), rB0
            #elif KB == 105
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4896(pA), rB1
            addpd rB1, rC21
            #if KB > 105
               movddup 2400(pB), rB1
            #elif KB == 105
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 105
               movaps 4912(pA), rA0
            #elif KB == 105
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 105
               movaps 4928(pA), rA1
            #elif KB == 105
               prefA(-128(pfA))
            #endif

            mulpd 4896(pA), rB2
            addpd rB2, rC22
            #if KB > 105
               movddup 2408(pB), rB2
            #elif KB == 105
               prefA((pfA))
            #endif

         #endif
         #if KB > 105
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4944(pA), rB0
            addpd rB0, rC20
            #if KB > 106
               movddup 2416(pB), rB0
            #elif KB == 106
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4944(pA), rB1
            addpd rB1, rC21
            #if KB > 106
               movddup 2424(pB), rB1
            #elif KB == 106
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 106
               movaps 4960(pA), rA0
            #elif KB == 106
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 106
               movaps 4976(pA), rA1
            #elif KB == 106
               prefA(-128(pfA))
            #endif

            mulpd 4944(pA), rB2
            addpd rB2, rC22
            #if KB > 106
               movddup 2432(pB), rB2
            #elif KB == 106
               prefA((pfA))
            #endif

         #endif
         #if KB > 106
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 4992(pA), rB0
            addpd rB0, rC20
            #if KB > 107
               movddup 2440(pB), rB0
            #elif KB == 107
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 4992(pA), rB1
            addpd rB1, rC21
            #if KB > 107
               movddup 2448(pB), rB1
            #elif KB == 107
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 107
               movaps 5008(pA), rA0
            #elif KB == 107
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 107
               movaps 5024(pA), rA1
            #elif KB == 107
               prefA(-128(pfA))
            #endif

            mulpd 4992(pA), rB2
            addpd rB2, rC22
            #if KB > 107
               movddup 2456(pB), rB2
            #elif KB == 107
               prefA((pfA))
            #endif

         #endif
         #if KB > 107
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5040(pA), rB0
            addpd rB0, rC20
            #if KB > 108
               movddup 2464(pB), rB0
            #elif KB == 108
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5040(pA), rB1
            addpd rB1, rC21
            #if KB > 108
               movddup 2472(pB), rB1
            #elif KB == 108
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 108
               movaps 5056(pA), rA0
            #elif KB == 108
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 108
               movaps 5072(pA), rA1
            #elif KB == 108
               prefA(-128(pfA))
            #endif

            mulpd 5040(pA), rB2
            addpd rB2, rC22
            #if KB > 108
               movddup 2480(pB), rB2
            #elif KB == 108
               prefA((pfA))
            #endif

         #endif
         #if KB > 108
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5088(pA), rB0
            addpd rB0, rC20
            #if KB > 109
               movddup 2488(pB), rB0
            #elif KB == 109
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5088(pA), rB1
            addpd rB1, rC21
            #if KB > 109
               movddup 2496(pB), rB1
            #elif KB == 109
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 109
               movaps 5104(pA), rA0
            #elif KB == 109
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 109
               movaps 5120(pA), rA1
            #elif KB == 109
               prefA(-128(pfA))
            #endif

            mulpd 5088(pA), rB2
            addpd rB2, rC22
            #if KB > 109
               movddup 2504(pB), rB2
            #elif KB == 109
               prefA((pfA))
            #endif

         #endif
         #if KB > 109
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5136(pA), rB0
            addpd rB0, rC20
            #if KB > 110
               movddup 2512(pB), rB0
            #elif KB == 110
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5136(pA), rB1
            addpd rB1, rC21
            #if KB > 110
               movddup 2520(pB), rB1
            #elif KB == 110
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 110
               movaps 5152(pA), rA0
            #elif KB == 110
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 110
               movaps 5168(pA), rA1
            #elif KB == 110
               prefA(-128(pfA))
            #endif

            mulpd 5136(pA), rB2
            addpd rB2, rC22
            #if KB > 110
               movddup 2528(pB), rB2
            #elif KB == 110
               prefA((pfA))
            #endif

         #endif
         #if KB > 110
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5184(pA), rB0
            addpd rB0, rC20
            #if KB > 111
               movddup 2536(pB), rB0
            #elif KB == 111
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5184(pA), rB1
            addpd rB1, rC21
            #if KB > 111
               movddup 2544(pB), rB1
            #elif KB == 111
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 111
               movaps 5200(pA), rA0
            #elif KB == 111
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 111
               movaps 5216(pA), rA1
            #elif KB == 111
               prefA(-128(pfA))
            #endif

            mulpd 5184(pA), rB2
            addpd rB2, rC22
            #if KB > 111
               movddup 2552(pB), rB2
            #elif KB == 111
               prefA((pfA))
            #endif

         #endif
         #if KB > 111
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5232(pA), rB0
            addpd rB0, rC20
            #if KB > 112
               movddup 2560(pB), rB0
            #elif KB == 112
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5232(pA), rB1
            addpd rB1, rC21
            #if KB > 112
               movddup 2568(pB), rB1
            #elif KB == 112
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 112
               movaps 5248(pA), rA0
            #elif KB == 112
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 112
               movaps 5264(pA), rA1
            #elif KB == 112
               prefA(-128(pfA))
            #endif

            mulpd 5232(pA), rB2
            addpd rB2, rC22
            #if KB > 112
               movddup 2576(pB), rB2
            #elif KB == 112
               prefA((pfA))
            #endif

         #endif
         #if KB > 112
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5280(pA), rB0
            addpd rB0, rC20
            #if KB > 113
               movddup 2584(pB), rB0
            #elif KB == 113
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5280(pA), rB1
            addpd rB1, rC21
            #if KB > 113
               movddup 2592(pB), rB1
            #elif KB == 113
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 113
               movaps 5296(pA), rA0
            #elif KB == 113
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 113
               movaps 5312(pA), rA1
            #elif KB == 113
               prefA(-128(pfA))
            #endif

            mulpd 5280(pA), rB2
            addpd rB2, rC22
            #if KB > 113
               movddup 2600(pB), rB2
            #elif KB == 113
               prefA((pfA))
            #endif

         #endif
         #if KB > 113
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5328(pA), rB0
            addpd rB0, rC20
            #if KB > 114
               movddup 2608(pB), rB0
            #elif KB == 114
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5328(pA), rB1
            addpd rB1, rC21
            #if KB > 114
               movddup 2616(pB), rB1
            #elif KB == 114
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 114
               movaps 5344(pA), rA0
            #elif KB == 114
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 114
               movaps 5360(pA), rA1
            #elif KB == 114
               prefA(-128(pfA))
            #endif

            mulpd 5328(pA), rB2
            addpd rB2, rC22
            #if KB > 114
               movddup 2624(pB), rB2
            #elif KB == 114
               prefA((pfA))
            #endif

         #endif
         #if KB > 114
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5376(pA), rB0
            addpd rB0, rC20
            #if KB > 115
               movddup 2632(pB), rB0
            #elif KB == 115
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5376(pA), rB1
            addpd rB1, rC21
            #if KB > 115
               movddup 2640(pB), rB1
            #elif KB == 115
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 115
               movaps 5392(pA), rA0
            #elif KB == 115
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 115
               movaps 5408(pA), rA1
            #elif KB == 115
               prefA(-128(pfA))
            #endif

            mulpd 5376(pA), rB2
            addpd rB2, rC22
            #if KB > 115
               movddup 2648(pB), rB2
            #elif KB == 115
               prefA((pfA))
            #endif

         #endif
         #if KB > 115
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5424(pA), rB0
            addpd rB0, rC20
            #if KB > 116
               movddup 2656(pB), rB0
            #elif KB == 116
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5424(pA), rB1
            addpd rB1, rC21
            #if KB > 116
               movddup 2664(pB), rB1
            #elif KB == 116
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 116
               movaps 5440(pA), rA0
            #elif KB == 116
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 116
               movaps 5456(pA), rA1
            #elif KB == 116
               prefA(-128(pfA))
            #endif

            mulpd 5424(pA), rB2
            addpd rB2, rC22
            #if KB > 116
               movddup 2672(pB), rB2
            #elif KB == 116
               prefA((pfA))
            #endif

         #endif
         #if KB > 116
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5472(pA), rB0
            addpd rB0, rC20
            #if KB > 117
               movddup 2680(pB), rB0
            #elif KB == 117
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5472(pA), rB1
            addpd rB1, rC21
            #if KB > 117
               movddup 2688(pB), rB1
            #elif KB == 117
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 117
               movaps 5488(pA), rA0
            #elif KB == 117
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 117
               movaps 5504(pA), rA1
            #elif KB == 117
               prefA(-128(pfA))
            #endif

            mulpd 5472(pA), rB2
            addpd rB2, rC22
            #if KB > 117
               movddup 2696(pB), rB2
            #elif KB == 117
               prefA((pfA))
            #endif

         #endif
         #if KB > 117
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5520(pA), rB0
            addpd rB0, rC20
            #if KB > 118
               movddup 2704(pB), rB0
            #elif KB == 118
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5520(pA), rB1
            addpd rB1, rC21
            #if KB > 118
               movddup 2712(pB), rB1
            #elif KB == 118
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 118
               movaps 5536(pA), rA0
            #elif KB == 118
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 118
               movaps 5552(pA), rA1
            #elif KB == 118
               prefA(-128(pfA))
            #endif

            mulpd 5520(pA), rB2
            addpd rB2, rC22
            #if KB > 118
               movddup 2720(pB), rB2
            #elif KB == 118
               prefA((pfA))
            #endif

         #endif
         #if KB > 118
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5568(pA), rB0
            addpd rB0, rC20
            #if KB > 119
               movddup 2728(pB), rB0
            #elif KB == 119
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5568(pA), rB1
            addpd rB1, rC21
            #if KB > 119
               movddup 2736(pB), rB1
            #elif KB == 119
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 119
               movaps 5584(pA), rA0
            #elif KB == 119
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 119
               movaps 5600(pA), rA1
            #elif KB == 119
               prefA(-128(pfA))
            #endif

            mulpd 5568(pA), rB2
            addpd rB2, rC22
            #if KB > 119
               movddup 2744(pB), rB2
            #elif KB == 119
               prefA((pfA))
            #endif

         #endif
         #if KB > 119
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5616(pA), rB0
            addpd rB0, rC20
            #if KB > 120
               movddup 2752(pB), rB0
            #elif KB == 120
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5616(pA), rB1
            addpd rB1, rC21
            #if KB > 120
               movddup 2760(pB), rB1
            #elif KB == 120
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 120
               movaps 5632(pA), rA0
            #elif KB == 120
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 120
               movaps 5648(pA), rA1
            #elif KB == 120
               prefA(-128(pfA))
            #endif

            mulpd 5616(pA), rB2
            addpd rB2, rC22
            #if KB > 120
               movddup 2768(pB), rB2
            #elif KB == 120
               prefA((pfA))
            #endif

         #endif
         #if KB > 120
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5664(pA), rB0
            addpd rB0, rC20
            #if KB > 121
               movddup 2776(pB), rB0
            #elif KB == 121
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5664(pA), rB1
            addpd rB1, rC21
            #if KB > 121
               movddup 2784(pB), rB1
            #elif KB == 121
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 121
               movaps 5680(pA), rA0
            #elif KB == 121
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 121
               movaps 5696(pA), rA1
            #elif KB == 121
               prefA(-128(pfA))
            #endif

            mulpd 5664(pA), rB2
            addpd rB2, rC22
            #if KB > 121
               movddup 2792(pB), rB2
            #elif KB == 121
               prefA((pfA))
            #endif

         #endif
         #if KB > 121
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5712(pA), rB0
            addpd rB0, rC20
            #if KB > 122
               movddup 2800(pB), rB0
            #elif KB == 122
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5712(pA), rB1
            addpd rB1, rC21
            #if KB > 122
               movddup 2808(pB), rB1
            #elif KB == 122
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 122
               movaps 5728(pA), rA0
            #elif KB == 122
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 122
               movaps 5744(pA), rA1
            #elif KB == 122
               prefA(-128(pfA))
            #endif

            mulpd 5712(pA), rB2
            addpd rB2, rC22
            #if KB > 122
               movddup 2816(pB), rB2
            #elif KB == 122
               prefA((pfA))
            #endif

         #endif
         #if KB > 122
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5760(pA), rB0
            addpd rB0, rC20
            #if KB > 123
               movddup 2824(pB), rB0
            #elif KB == 123
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5760(pA), rB1
            addpd rB1, rC21
            #if KB > 123
               movddup 2832(pB), rB1
            #elif KB == 123
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 123
               movaps 5776(pA), rA0
            #elif KB == 123
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 123
               movaps 5792(pA), rA1
            #elif KB == 123
               prefA(-128(pfA))
            #endif

            mulpd 5760(pA), rB2
            addpd rB2, rC22
            #if KB > 123
               movddup 2840(pB), rB2
            #elif KB == 123
               prefA((pfA))
            #endif

         #endif
         #if KB > 123
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5808(pA), rB0
            addpd rB0, rC20
            #if KB > 124
               movddup 2848(pB), rB0
            #elif KB == 124
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5808(pA), rB1
            addpd rB1, rC21
            #if KB > 124
               movddup 2856(pB), rB1
            #elif KB == 124
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 124
               movaps 5824(pA), rA0
            #elif KB == 124
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 124
               movaps 5840(pA), rA1
            #elif KB == 124
               prefA(-128(pfA))
            #endif

            mulpd 5808(pA), rB2
            addpd rB2, rC22
            #if KB > 124
               movddup 2864(pB), rB2
            #elif KB == 124
               prefA((pfA))
            #endif

         #endif
         #if KB > 124
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5856(pA), rB0
            addpd rB0, rC20
            #if KB > 125
               movddup 2872(pB), rB0
            #elif KB == 125
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5856(pA), rB1
            addpd rB1, rC21
            #if KB > 125
               movddup 2880(pB), rB1
            #elif KB == 125
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 125
               movaps 5872(pA), rA0
            #elif KB == 125
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 125
               movaps 5888(pA), rA1
            #elif KB == 125
               prefA(-128(pfA))
            #endif

            mulpd 5856(pA), rB2
            addpd rB2, rC22
            #if KB > 125
               movddup 2888(pB), rB2
            #elif KB == 125
               prefA((pfA))
            #endif

         #endif
         #if KB > 125
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5904(pA), rB0
            addpd rB0, rC20
            #if KB > 126
               movddup 2896(pB), rB0
            #elif KB == 126
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5904(pA), rB1
            addpd rB1, rC21
            #if KB > 126
               movddup 2904(pB), rB1
            #elif KB == 126
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 126
               movaps 5920(pA), rA0
            #elif KB == 126
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 126
               movaps 5936(pA), rA1
            #elif KB == 126
               prefA(-128(pfA))
            #endif

            mulpd 5904(pA), rB2
            addpd rB2, rC22
            #if KB > 126
               movddup 2912(pB), rB2
            #elif KB == 126
               prefA((pfA))
            #endif

         #endif
         #if KB > 126
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 5952(pA), rB0
            addpd rB0, rC20
            #if KB > 127
               movddup 2920(pB), rB0
            #elif KB == 127
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 5952(pA), rB1
            addpd rB1, rC21
            #if KB > 127
               movddup 2928(pB), rB1
            #elif KB == 127
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 127
               movaps 5968(pA), rA0
            #elif KB == 127
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 127
               movaps 5984(pA), rA1
            #elif KB == 127
               prefA(-128(pfA))
            #endif

            mulpd 5952(pA), rB2
            addpd rB2, rC22
            #if KB > 127
               movddup 2936(pB), rB2
            #elif KB == 127
               prefA((pfA))
            #endif

         #endif
         #if KB > 127
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6000(pA), rB0
            addpd rB0, rC20
            #if KB > 128
               movddup 2944(pB), rB0
            #elif KB == 128
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6000(pA), rB1
            addpd rB1, rC21
            #if KB > 128
               movddup 2952(pB), rB1
            #elif KB == 128
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 128
               movaps 6016(pA), rA0
            #elif KB == 128
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 128
               movaps 6032(pA), rA1
            #elif KB == 128
               prefA(-128(pfA))
            #endif

            mulpd 6000(pA), rB2
            addpd rB2, rC22
            #if KB > 128
               movddup 2960(pB), rB2
            #elif KB == 128
               prefA((pfA))
            #endif

         #endif
         #if KB > 128
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6048(pA), rB0
            addpd rB0, rC20
            #if KB > 129
               movddup 2968(pB), rB0
            #elif KB == 129
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6048(pA), rB1
            addpd rB1, rC21
            #if KB > 129
               movddup 2976(pB), rB1
            #elif KB == 129
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 129
               movaps 6064(pA), rA0
            #elif KB == 129
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 129
               movaps 6080(pA), rA1
            #elif KB == 129
               prefA(-128(pfA))
            #endif

            mulpd 6048(pA), rB2
            addpd rB2, rC22
            #if KB > 129
               movddup 2984(pB), rB2
            #elif KB == 129
               prefA((pfA))
            #endif

         #endif
         #if KB > 129
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6096(pA), rB0
            addpd rB0, rC20
            #if KB > 130
               movddup 2992(pB), rB0
            #elif KB == 130
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6096(pA), rB1
            addpd rB1, rC21
            #if KB > 130
               movddup 3000(pB), rB1
            #elif KB == 130
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 130
               movaps 6112(pA), rA0
            #elif KB == 130
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 130
               movaps 6128(pA), rA1
            #elif KB == 130
               prefA(-128(pfA))
            #endif

            mulpd 6096(pA), rB2
            addpd rB2, rC22
            #if KB > 130
               movddup 3008(pB), rB2
            #elif KB == 130
               prefA((pfA))
            #endif

         #endif
         #if KB > 130
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6144(pA), rB0
            addpd rB0, rC20
            #if KB > 131
               movddup 3016(pB), rB0
            #elif KB == 131
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6144(pA), rB1
            addpd rB1, rC21
            #if KB > 131
               movddup 3024(pB), rB1
            #elif KB == 131
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 131
               movaps 6160(pA), rA0
            #elif KB == 131
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 131
               movaps 6176(pA), rA1
            #elif KB == 131
               prefA(-128(pfA))
            #endif

            mulpd 6144(pA), rB2
            addpd rB2, rC22
            #if KB > 131
               movddup 3032(pB), rB2
            #elif KB == 131
               prefA((pfA))
            #endif

         #endif
         #if KB > 131
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6192(pA), rB0
            addpd rB0, rC20
            #if KB > 132
               movddup 3040(pB), rB0
            #elif KB == 132
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6192(pA), rB1
            addpd rB1, rC21
            #if KB > 132
               movddup 3048(pB), rB1
            #elif KB == 132
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 132
               movaps 6208(pA), rA0
            #elif KB == 132
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 132
               movaps 6224(pA), rA1
            #elif KB == 132
               prefA(-128(pfA))
            #endif

            mulpd 6192(pA), rB2
            addpd rB2, rC22
            #if KB > 132
               movddup 3056(pB), rB2
            #elif KB == 132
               prefA((pfA))
            #endif

         #endif
         #if KB > 132
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6240(pA), rB0
            addpd rB0, rC20
            #if KB > 133
               movddup 3064(pB), rB0
            #elif KB == 133
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6240(pA), rB1
            addpd rB1, rC21
            #if KB > 133
               movddup 3072(pB), rB1
            #elif KB == 133
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 133
               movaps 6256(pA), rA0
            #elif KB == 133
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 133
               movaps 6272(pA), rA1
            #elif KB == 133
               prefA(-128(pfA))
            #endif

            mulpd 6240(pA), rB2
            addpd rB2, rC22
            #if KB > 133
               movddup 3080(pB), rB2
            #elif KB == 133
               prefA((pfA))
            #endif

         #endif
         #if KB > 133
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6288(pA), rB0
            addpd rB0, rC20
            #if KB > 134
               movddup 3088(pB), rB0
            #elif KB == 134
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6288(pA), rB1
            addpd rB1, rC21
            #if KB > 134
               movddup 3096(pB), rB1
            #elif KB == 134
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 134
               movaps 6304(pA), rA0
            #elif KB == 134
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 134
               movaps 6320(pA), rA1
            #elif KB == 134
               prefA(-128(pfA))
            #endif

            mulpd 6288(pA), rB2
            addpd rB2, rC22
            #if KB > 134
               movddup 3104(pB), rB2
            #elif KB == 134
               prefA((pfA))
            #endif

         #endif
         #if KB > 134
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6336(pA), rB0
            addpd rB0, rC20
            #if KB > 135
               movddup 3112(pB), rB0
            #elif KB == 135
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6336(pA), rB1
            addpd rB1, rC21
            #if KB > 135
               movddup 3120(pB), rB1
            #elif KB == 135
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 135
               movaps 6352(pA), rA0
            #elif KB == 135
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 135
               movaps 6368(pA), rA1
            #elif KB == 135
               prefA(-128(pfA))
            #endif

            mulpd 6336(pA), rB2
            addpd rB2, rC22
            #if KB > 135
               movddup 3128(pB), rB2
            #elif KB == 135
               prefA((pfA))
            #endif

         #endif
         #if KB > 135
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6384(pA), rB0
            addpd rB0, rC20
            #if KB > 136
               movddup 3136(pB), rB0
            #elif KB == 136
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6384(pA), rB1
            addpd rB1, rC21
            #if KB > 136
               movddup 3144(pB), rB1
            #elif KB == 136
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 136
               movaps 6400(pA), rA0
            #elif KB == 136
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 136
               movaps 6416(pA), rA1
            #elif KB == 136
               prefA(-128(pfA))
            #endif

            mulpd 6384(pA), rB2
            addpd rB2, rC22
            #if KB > 136
               movddup 3152(pB), rB2
            #elif KB == 136
               prefA((pfA))
            #endif

         #endif
         #if KB > 136
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6432(pA), rB0
            addpd rB0, rC20
            #if KB > 137
               movddup 3160(pB), rB0
            #elif KB == 137
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6432(pA), rB1
            addpd rB1, rC21
            #if KB > 137
               movddup 3168(pB), rB1
            #elif KB == 137
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 137
               movaps 6448(pA), rA0
            #elif KB == 137
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 137
               movaps 6464(pA), rA1
            #elif KB == 137
               prefA(-128(pfA))
            #endif

            mulpd 6432(pA), rB2
            addpd rB2, rC22
            #if KB > 137
               movddup 3176(pB), rB2
            #elif KB == 137
               prefA((pfA))
            #endif

         #endif
         #if KB > 137
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6480(pA), rB0
            addpd rB0, rC20
            #if KB > 138
               movddup 3184(pB), rB0
            #elif KB == 138
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6480(pA), rB1
            addpd rB1, rC21
            #if KB > 138
               movddup 3192(pB), rB1
            #elif KB == 138
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 138
               movaps 6496(pA), rA0
            #elif KB == 138
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 138
               movaps 6512(pA), rA1
            #elif KB == 138
               prefA(-128(pfA))
            #endif

            mulpd 6480(pA), rB2
            addpd rB2, rC22
            #if KB > 138
               movddup 3200(pB), rB2
            #elif KB == 138
               prefA((pfA))
            #endif

         #endif
         #if KB > 138
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6528(pA), rB0
            addpd rB0, rC20
            #if KB > 139
               movddup 3208(pB), rB0
            #elif KB == 139
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6528(pA), rB1
            addpd rB1, rC21
            #if KB > 139
               movddup 3216(pB), rB1
            #elif KB == 139
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 139
               movaps 6544(pA), rA0
            #elif KB == 139
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 139
               movaps 6560(pA), rA1
            #elif KB == 139
               prefA(-128(pfA))
            #endif

            mulpd 6528(pA), rB2
            addpd rB2, rC22
            #if KB > 139
               movddup 3224(pB), rB2
            #elif KB == 139
               prefA((pfA))
            #endif

         #endif
         #if KB > 139
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6576(pA), rB0
            addpd rB0, rC20
            #if KB > 140
               movddup 3232(pB), rB0
            #elif KB == 140
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6576(pA), rB1
            addpd rB1, rC21
            #if KB > 140
               movddup 3240(pB), rB1
            #elif KB == 140
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 140
               movaps 6592(pA), rA0
            #elif KB == 140
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 140
               movaps 6608(pA), rA1
            #elif KB == 140
               prefA(-128(pfA))
            #endif

            mulpd 6576(pA), rB2
            addpd rB2, rC22
            #if KB > 140
               movddup 3248(pB), rB2
            #elif KB == 140
               prefA((pfA))
            #endif

         #endif
         #if KB > 140
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6624(pA), rB0
            addpd rB0, rC20
            #if KB > 141
               movddup 3256(pB), rB0
            #elif KB == 141
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6624(pA), rB1
            addpd rB1, rC21
            #if KB > 141
               movddup 3264(pB), rB1
            #elif KB == 141
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 141
               movaps 6640(pA), rA0
            #elif KB == 141
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 141
               movaps 6656(pA), rA1
            #elif KB == 141
               prefA(-128(pfA))
            #endif

            mulpd 6624(pA), rB2
            addpd rB2, rC22
            #if KB > 141
               movddup 3272(pB), rB2
            #elif KB == 141
               prefA((pfA))
            #endif

         #endif
         #if KB > 141
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6672(pA), rB0
            addpd rB0, rC20
            #if KB > 142
               movddup 3280(pB), rB0
            #elif KB == 142
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6672(pA), rB1
            addpd rB1, rC21
            #if KB > 142
               movddup 3288(pB), rB1
            #elif KB == 142
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 142
               movaps 6688(pA), rA0
            #elif KB == 142
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 142
               movaps 6704(pA), rA1
            #elif KB == 142
               prefA(-128(pfA))
            #endif

            mulpd 6672(pA), rB2
            addpd rB2, rC22
            #if KB > 142
               movddup 3296(pB), rB2
            #elif KB == 142
               prefA((pfA))
            #endif

         #endif
         #if KB > 142
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6720(pA), rB0
            addpd rB0, rC20
            #if KB > 143
               movddup 3304(pB), rB0
            #elif KB == 143
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6720(pA), rB1
            addpd rB1, rC21
            #if KB > 143
               movddup 3312(pB), rB1
            #elif KB == 143
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 143
               movaps 6736(pA), rA0
            #elif KB == 143
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 143
               movaps 6752(pA), rA1
            #elif KB == 143
               prefA(-128(pfA))
            #endif

            mulpd 6720(pA), rB2
            addpd rB2, rC22
            #if KB > 143
               movddup 3320(pB), rB2
            #elif KB == 143
               prefA((pfA))
            #endif

         #endif
         #if KB > 143
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6768(pA), rB0
            addpd rB0, rC20
            #if KB > 144
               movddup 3328(pB), rB0
            #elif KB == 144
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6768(pA), rB1
            addpd rB1, rC21
            #if KB > 144
               movddup 3336(pB), rB1
            #elif KB == 144
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 144
               movaps 6784(pA), rA0
            #elif KB == 144
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 144
               movaps 6800(pA), rA1
            #elif KB == 144
               prefA(-128(pfA))
            #endif

            mulpd 6768(pA), rB2
            addpd rB2, rC22
            #if KB > 144
               movddup 3344(pB), rB2
            #elif KB == 144
               prefA((pfA))
            #endif

         #endif
         #if KB > 144
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6816(pA), rB0
            addpd rB0, rC20
            #if KB > 145
               movddup 3352(pB), rB0
            #elif KB == 145
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6816(pA), rB1
            addpd rB1, rC21
            #if KB > 145
               movddup 3360(pB), rB1
            #elif KB == 145
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 145
               movaps 6832(pA), rA0
            #elif KB == 145
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 145
               movaps 6848(pA), rA1
            #elif KB == 145
               prefA(-128(pfA))
            #endif

            mulpd 6816(pA), rB2
            addpd rB2, rC22
            #if KB > 145
               movddup 3368(pB), rB2
            #elif KB == 145
               prefA((pfA))
            #endif

         #endif
         #if KB > 145
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6864(pA), rB0
            addpd rB0, rC20
            #if KB > 146
               movddup 3376(pB), rB0
            #elif KB == 146
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6864(pA), rB1
            addpd rB1, rC21
            #if KB > 146
               movddup 3384(pB), rB1
            #elif KB == 146
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 146
               movaps 6880(pA), rA0
            #elif KB == 146
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 146
               movaps 6896(pA), rA1
            #elif KB == 146
               prefA(-128(pfA))
            #endif

            mulpd 6864(pA), rB2
            addpd rB2, rC22
            #if KB > 146
               movddup 3392(pB), rB2
            #elif KB == 146
               prefA((pfA))
            #endif

         #endif
         #if KB > 146
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6912(pA), rB0
            addpd rB0, rC20
            #if KB > 147
               movddup 3400(pB), rB0
            #elif KB == 147
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6912(pA), rB1
            addpd rB1, rC21
            #if KB > 147
               movddup 3408(pB), rB1
            #elif KB == 147
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 147
               movaps 6928(pA), rA0
            #elif KB == 147
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 147
               movaps 6944(pA), rA1
            #elif KB == 147
               prefA(-128(pfA))
            #endif

            mulpd 6912(pA), rB2
            addpd rB2, rC22
            #if KB > 147
               movddup 3416(pB), rB2
            #elif KB == 147
               prefA((pfA))
            #endif

         #endif
         #if KB > 147
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 6960(pA), rB0
            addpd rB0, rC20
            #if KB > 148
               movddup 3424(pB), rB0
            #elif KB == 148
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 6960(pA), rB1
            addpd rB1, rC21
            #if KB > 148
               movddup 3432(pB), rB1
            #elif KB == 148
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 148
               movaps 6976(pA), rA0
            #elif KB == 148
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 148
               movaps 6992(pA), rA1
            #elif KB == 148
               prefA(-128(pfA))
            #endif

            mulpd 6960(pA), rB2
            addpd rB2, rC22
            #if KB > 148
               movddup 3440(pB), rB2
            #elif KB == 148
               prefA((pfA))
            #endif

         #endif
         #if KB > 148
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7008(pA), rB0
            addpd rB0, rC20
            #if KB > 149
               movddup 3448(pB), rB0
            #elif KB == 149
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7008(pA), rB1
            addpd rB1, rC21
            #if KB > 149
               movddup 3456(pB), rB1
            #elif KB == 149
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 149
               movaps 7024(pA), rA0
            #elif KB == 149
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 149
               movaps 7040(pA), rA1
            #elif KB == 149
               prefA(-128(pfA))
            #endif

            mulpd 7008(pA), rB2
            addpd rB2, rC22
            #if KB > 149
               movddup 3464(pB), rB2
            #elif KB == 149
               prefA((pfA))
            #endif

         #endif
         #if KB > 149
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7056(pA), rB0
            addpd rB0, rC20
            #if KB > 150
               movddup 3472(pB), rB0
            #elif KB == 150
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7056(pA), rB1
            addpd rB1, rC21
            #if KB > 150
               movddup 3480(pB), rB1
            #elif KB == 150
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 150
               movaps 7072(pA), rA0
            #elif KB == 150
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 150
               movaps 7088(pA), rA1
            #elif KB == 150
               prefA(-128(pfA))
            #endif

            mulpd 7056(pA), rB2
            addpd rB2, rC22
            #if KB > 150
               movddup 3488(pB), rB2
            #elif KB == 150
               prefA((pfA))
            #endif

         #endif
         #if KB > 150
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7104(pA), rB0
            addpd rB0, rC20
            #if KB > 151
               movddup 3496(pB), rB0
            #elif KB == 151
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7104(pA), rB1
            addpd rB1, rC21
            #if KB > 151
               movddup 3504(pB), rB1
            #elif KB == 151
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 151
               movaps 7120(pA), rA0
            #elif KB == 151
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 151
               movaps 7136(pA), rA1
            #elif KB == 151
               prefA(-128(pfA))
            #endif

            mulpd 7104(pA), rB2
            addpd rB2, rC22
            #if KB > 151
               movddup 3512(pB), rB2
            #elif KB == 151
               prefA((pfA))
            #endif

         #endif
         #if KB > 151
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7152(pA), rB0
            addpd rB0, rC20
            #if KB > 152
               movddup 3520(pB), rB0
            #elif KB == 152
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7152(pA), rB1
            addpd rB1, rC21
            #if KB > 152
               movddup 3528(pB), rB1
            #elif KB == 152
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 152
               movaps 7168(pA), rA0
            #elif KB == 152
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 152
               movaps 7184(pA), rA1
            #elif KB == 152
               prefA(-128(pfA))
            #endif

            mulpd 7152(pA), rB2
            addpd rB2, rC22
            #if KB > 152
               movddup 3536(pB), rB2
            #elif KB == 152
               prefA((pfA))
            #endif

         #endif
         #if KB > 152
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7200(pA), rB0
            addpd rB0, rC20
            #if KB > 153
               movddup 3544(pB), rB0
            #elif KB == 153
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7200(pA), rB1
            addpd rB1, rC21
            #if KB > 153
               movddup 3552(pB), rB1
            #elif KB == 153
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 153
               movaps 7216(pA), rA0
            #elif KB == 153
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 153
               movaps 7232(pA), rA1
            #elif KB == 153
               prefA(-128(pfA))
            #endif

            mulpd 7200(pA), rB2
            addpd rB2, rC22
            #if KB > 153
               movddup 3560(pB), rB2
            #elif KB == 153
               prefA((pfA))
            #endif

         #endif
         #if KB > 153
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7248(pA), rB0
            addpd rB0, rC20
            #if KB > 154
               movddup 3568(pB), rB0
            #elif KB == 154
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7248(pA), rB1
            addpd rB1, rC21
            #if KB > 154
               movddup 3576(pB), rB1
            #elif KB == 154
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 154
               movaps 7264(pA), rA0
            #elif KB == 154
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 154
               movaps 7280(pA), rA1
            #elif KB == 154
               prefA(-128(pfA))
            #endif

            mulpd 7248(pA), rB2
            addpd rB2, rC22
            #if KB > 154
               movddup 3584(pB), rB2
            #elif KB == 154
               prefA((pfA))
            #endif

         #endif
         #if KB > 154
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7296(pA), rB0
            addpd rB0, rC20
            #if KB > 155
               movddup 3592(pB), rB0
            #elif KB == 155
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7296(pA), rB1
            addpd rB1, rC21
            #if KB > 155
               movddup 3600(pB), rB1
            #elif KB == 155
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 155
               movaps 7312(pA), rA0
            #elif KB == 155
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 155
               movaps 7328(pA), rA1
            #elif KB == 155
               prefA(-128(pfA))
            #endif

            mulpd 7296(pA), rB2
            addpd rB2, rC22
            #if KB > 155
               movddup 3608(pB), rB2
            #elif KB == 155
               prefA((pfA))
            #endif

         #endif
         #if KB > 155
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7344(pA), rB0
            addpd rB0, rC20
            #if KB > 156
               movddup 3616(pB), rB0
            #elif KB == 156
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7344(pA), rB1
            addpd rB1, rC21
            #if KB > 156
               movddup 3624(pB), rB1
            #elif KB == 156
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 156
               movaps 7360(pA), rA0
            #elif KB == 156
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 156
               movaps 7376(pA), rA1
            #elif KB == 156
               prefA(-128(pfA))
            #endif

            mulpd 7344(pA), rB2
            addpd rB2, rC22
            #if KB > 156
               movddup 3632(pB), rB2
            #elif KB == 156
               prefA((pfA))
            #endif

         #endif
         #if KB > 156
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7392(pA), rB0
            addpd rB0, rC20
            #if KB > 157
               movddup 3640(pB), rB0
            #elif KB == 157
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7392(pA), rB1
            addpd rB1, rC21
            #if KB > 157
               movddup 3648(pB), rB1
            #elif KB == 157
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 157
               movaps 7408(pA), rA0
            #elif KB == 157
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 157
               movaps 7424(pA), rA1
            #elif KB == 157
               prefA(-128(pfA))
            #endif

            mulpd 7392(pA), rB2
            addpd rB2, rC22
            #if KB > 157
               movddup 3656(pB), rB2
            #elif KB == 157
               prefA((pfA))
            #endif

         #endif
         #if KB > 157
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7440(pA), rB0
            addpd rB0, rC20
            #if KB > 158
               movddup 3664(pB), rB0
            #elif KB == 158
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7440(pA), rB1
            addpd rB1, rC21
            #if KB > 158
               movddup 3672(pB), rB1
            #elif KB == 158
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 158
               movaps 7456(pA), rA0
            #elif KB == 158
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 158
               movaps 7472(pA), rA1
            #elif KB == 158
               prefA(-128(pfA))
            #endif

            mulpd 7440(pA), rB2
            addpd rB2, rC22
            #if KB > 158
               movddup 3680(pB), rB2
            #elif KB == 158
               prefA((pfA))
            #endif

         #endif
         #if KB > 158
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7488(pA), rB0
            addpd rB0, rC20
            #if KB > 159
               movddup 3688(pB), rB0
            #elif KB == 159
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7488(pA), rB1
            addpd rB1, rC21
            #if KB > 159
               movddup 3696(pB), rB1
            #elif KB == 159
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 159
               movaps 7504(pA), rA0
            #elif KB == 159
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 159
               movaps 7520(pA), rA1
            #elif KB == 159
               prefA(-128(pfA))
            #endif

            mulpd 7488(pA), rB2
            addpd rB2, rC22
            #if KB > 159
               movddup 3704(pB), rB2
            #elif KB == 159
               prefA((pfA))
            #endif

         #endif
         #if KB > 159
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7536(pA), rB0
            addpd rB0, rC20
            #if KB > 160
               movddup 3712(pB), rB0
            #elif KB == 160
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7536(pA), rB1
            addpd rB1, rC21
            #if KB > 160
               movddup 3720(pB), rB1
            #elif KB == 160
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 160
               movaps 7552(pA), rA0
            #elif KB == 160
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 160
               movaps 7568(pA), rA1
            #elif KB == 160
               prefA(-128(pfA))
            #endif

            mulpd 7536(pA), rB2
            addpd rB2, rC22
            #if KB > 160
               movddup 3728(pB), rB2
            #elif KB == 160
               prefA((pfA))
            #endif

         #endif
         #if KB > 160
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7584(pA), rB0
            addpd rB0, rC20
            #if KB > 161
               movddup 3736(pB), rB0
            #elif KB == 161
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7584(pA), rB1
            addpd rB1, rC21
            #if KB > 161
               movddup 3744(pB), rB1
            #elif KB == 161
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 161
               movaps 7600(pA), rA0
            #elif KB == 161
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 161
               movaps 7616(pA), rA1
            #elif KB == 161
               prefA(-128(pfA))
            #endif

            mulpd 7584(pA), rB2
            addpd rB2, rC22
            #if KB > 161
               movddup 3752(pB), rB2
            #elif KB == 161
               prefA((pfA))
            #endif

         #endif
         #if KB > 161
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7632(pA), rB0
            addpd rB0, rC20
            #if KB > 162
               movddup 3760(pB), rB0
            #elif KB == 162
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7632(pA), rB1
            addpd rB1, rC21
            #if KB > 162
               movddup 3768(pB), rB1
            #elif KB == 162
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 162
               movaps 7648(pA), rA0
            #elif KB == 162
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 162
               movaps 7664(pA), rA1
            #elif KB == 162
               prefA(-128(pfA))
            #endif

            mulpd 7632(pA), rB2
            addpd rB2, rC22
            #if KB > 162
               movddup 3776(pB), rB2
            #elif KB == 162
               prefA((pfA))
            #endif

         #endif
         #if KB > 162
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7680(pA), rB0
            addpd rB0, rC20
            #if KB > 163
               movddup 3784(pB), rB0
            #elif KB == 163
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7680(pA), rB1
            addpd rB1, rC21
            #if KB > 163
               movddup 3792(pB), rB1
            #elif KB == 163
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 163
               movaps 7696(pA), rA0
            #elif KB == 163
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 163
               movaps 7712(pA), rA1
            #elif KB == 163
               prefA(-128(pfA))
            #endif

            mulpd 7680(pA), rB2
            addpd rB2, rC22
            #if KB > 163
               movddup 3800(pB), rB2
            #elif KB == 163
               prefA((pfA))
            #endif

         #endif
         #if KB > 163
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7728(pA), rB0
            addpd rB0, rC20
            #if KB > 164
               movddup 3808(pB), rB0
            #elif KB == 164
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7728(pA), rB1
            addpd rB1, rC21
            #if KB > 164
               movddup 3816(pB), rB1
            #elif KB == 164
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 164
               movaps 7744(pA), rA0
            #elif KB == 164
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 164
               movaps 7760(pA), rA1
            #elif KB == 164
               prefA(-128(pfA))
            #endif

            mulpd 7728(pA), rB2
            addpd rB2, rC22
            #if KB > 164
               movddup 3824(pB), rB2
            #elif KB == 164
               prefA((pfA))
            #endif

         #endif
         #if KB > 164
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7776(pA), rB0
            addpd rB0, rC20
            #if KB > 165
               movddup 3832(pB), rB0
            #elif KB == 165
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7776(pA), rB1
            addpd rB1, rC21
            #if KB > 165
               movddup 3840(pB), rB1
            #elif KB == 165
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 165
               movaps 7792(pA), rA0
            #elif KB == 165
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 165
               movaps 7808(pA), rA1
            #elif KB == 165
               prefA(-128(pfA))
            #endif

            mulpd 7776(pA), rB2
            addpd rB2, rC22
            #if KB > 165
               movddup 3848(pB), rB2
            #elif KB == 165
               prefA((pfA))
            #endif

         #endif
         #if KB > 165
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7824(pA), rB0
            addpd rB0, rC20
            #if KB > 166
               movddup 3856(pB), rB0
            #elif KB == 166
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7824(pA), rB1
            addpd rB1, rC21
            #if KB > 166
               movddup 3864(pB), rB1
            #elif KB == 166
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 166
               movaps 7840(pA), rA0
            #elif KB == 166
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 166
               movaps 7856(pA), rA1
            #elif KB == 166
               prefA(-128(pfA))
            #endif

            mulpd 7824(pA), rB2
            addpd rB2, rC22
            #if KB > 166
               movddup 3872(pB), rB2
            #elif KB == 166
               prefA((pfA))
            #endif

         #endif
         #if KB > 166
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7872(pA), rB0
            addpd rB0, rC20
            #if KB > 167
               movddup 3880(pB), rB0
            #elif KB == 167
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7872(pA), rB1
            addpd rB1, rC21
            #if KB > 167
               movddup 3888(pB), rB1
            #elif KB == 167
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 167
               movaps 7888(pA), rA0
            #elif KB == 167
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 167
               movaps 7904(pA), rA1
            #elif KB == 167
               prefA(-128(pfA))
            #endif

            mulpd 7872(pA), rB2
            addpd rB2, rC22
            #if KB > 167
               movddup 3896(pB), rB2
            #elif KB == 167
               prefA((pfA))
            #endif

         #endif
         #if KB > 167
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7920(pA), rB0
            addpd rB0, rC20
            #if KB > 168
               movddup 3904(pB), rB0
            #elif KB == 168
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7920(pA), rB1
            addpd rB1, rC21
            #if KB > 168
               movddup 3912(pB), rB1
            #elif KB == 168
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 168
               movaps 7936(pA), rA0
            #elif KB == 168
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 168
               movaps 7952(pA), rA1
            #elif KB == 168
               prefA(-128(pfA))
            #endif

            mulpd 7920(pA), rB2
            addpd rB2, rC22
            #if KB > 168
               movddup 3920(pB), rB2
            #elif KB == 168
               prefA((pfA))
            #endif

         #endif
         #if KB > 168
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 7968(pA), rB0
            addpd rB0, rC20
            #if KB > 169
               movddup 3928(pB), rB0
            #elif KB == 169
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 7968(pA), rB1
            addpd rB1, rC21
            #if KB > 169
               movddup 3936(pB), rB1
            #elif KB == 169
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 169
               movaps 7984(pA), rA0
            #elif KB == 169
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 169
               movaps 8000(pA), rA1
            #elif KB == 169
               prefA(-128(pfA))
            #endif

            mulpd 7968(pA), rB2
            addpd rB2, rC22
            #if KB > 169
               movddup 3944(pB), rB2
            #elif KB == 169
               prefA((pfA))
            #endif

         #endif
         #if KB > 169
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8016(pA), rB0
            addpd rB0, rC20
            #if KB > 170
               movddup 3952(pB), rB0
            #elif KB == 170
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8016(pA), rB1
            addpd rB1, rC21
            #if KB > 170
               movddup 3960(pB), rB1
            #elif KB == 170
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 170
               movaps 8032(pA), rA0
            #elif KB == 170
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 170
               movaps 8048(pA), rA1
            #elif KB == 170
               prefA(-128(pfA))
            #endif

            mulpd 8016(pA), rB2
            addpd rB2, rC22
            #if KB > 170
               movddup 3968(pB), rB2
            #elif KB == 170
               prefA((pfA))
            #endif

         #endif
         #if KB > 170
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8064(pA), rB0
            addpd rB0, rC20
            #if KB > 171
               movddup 3976(pB), rB0
            #elif KB == 171
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8064(pA), rB1
            addpd rB1, rC21
            #if KB > 171
               movddup 3984(pB), rB1
            #elif KB == 171
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 171
               movaps 8080(pA), rA0
            #elif KB == 171
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 171
               movaps 8096(pA), rA1
            #elif KB == 171
               prefA(-128(pfA))
            #endif

            mulpd 8064(pA), rB2
            addpd rB2, rC22
            #if KB > 171
               movddup 3992(pB), rB2
            #elif KB == 171
               prefA((pfA))
            #endif

         #endif
         #if KB > 171
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8112(pA), rB0
            addpd rB0, rC20
            #if KB > 172
               movddup 4000(pB), rB0
            #elif KB == 172
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8112(pA), rB1
            addpd rB1, rC21
            #if KB > 172
               movddup 4008(pB), rB1
            #elif KB == 172
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 172
               movaps 8128(pA), rA0
            #elif KB == 172
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 172
               movaps 8144(pA), rA1
            #elif KB == 172
               prefA(-128(pfA))
            #endif

            mulpd 8112(pA), rB2
            addpd rB2, rC22
            #if KB > 172
               movddup 4016(pB), rB2
            #elif KB == 172
               prefA((pfA))
            #endif

         #endif
         #if KB > 172
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8160(pA), rB0
            addpd rB0, rC20
            #if KB > 173
               movddup 4024(pB), rB0
            #elif KB == 173
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8160(pA), rB1
            addpd rB1, rC21
            #if KB > 173
               movddup 4032(pB), rB1
            #elif KB == 173
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 173
               movaps 8176(pA), rA0
            #elif KB == 173
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 173
               movaps 8192(pA), rA1
            #elif KB == 173
               prefA(-128(pfA))
            #endif

            mulpd 8160(pA), rB2
            addpd rB2, rC22
            #if KB > 173
               movddup 4040(pB), rB2
            #elif KB == 173
               prefA((pfA))
            #endif

         #endif
         #if KB > 173
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8208(pA), rB0
            addpd rB0, rC20
            #if KB > 174
               movddup 4048(pB), rB0
            #elif KB == 174
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8208(pA), rB1
            addpd rB1, rC21
            #if KB > 174
               movddup 4056(pB), rB1
            #elif KB == 174
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 174
               movaps 8224(pA), rA0
            #elif KB == 174
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 174
               movaps 8240(pA), rA1
            #elif KB == 174
               prefA(-128(pfA))
            #endif

            mulpd 8208(pA), rB2
            addpd rB2, rC22
            #if KB > 174
               movddup 4064(pB), rB2
            #elif KB == 174
               prefA((pfA))
            #endif

         #endif
         #if KB > 174
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8256(pA), rB0
            addpd rB0, rC20
            #if KB > 175
               movddup 4072(pB), rB0
            #elif KB == 175
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8256(pA), rB1
            addpd rB1, rC21
            #if KB > 175
               movddup 4080(pB), rB1
            #elif KB == 175
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 175
               movaps 8272(pA), rA0
            #elif KB == 175
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 175
               movaps 8288(pA), rA1
            #elif KB == 175
               prefA(-128(pfA))
            #endif

            mulpd 8256(pA), rB2
            addpd rB2, rC22
            #if KB > 175
               movddup 4088(pB), rB2
            #elif KB == 175
               prefA((pfA))
            #endif

         #endif
         #if KB > 175
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8304(pA), rB0
            addpd rB0, rC20
            #if KB > 176
               movddup 4096(pB), rB0
            #elif KB == 176
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8304(pA), rB1
            addpd rB1, rC21
            #if KB > 176
               movddup 4104(pB), rB1
            #elif KB == 176
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 176
               movaps 8320(pA), rA0
            #elif KB == 176
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 176
               movaps 8336(pA), rA1
            #elif KB == 176
               prefA(-128(pfA))
            #endif

            mulpd 8304(pA), rB2
            addpd rB2, rC22
            #if KB > 176
               movddup 4112(pB), rB2
            #elif KB == 176
               prefA((pfA))
            #endif

         #endif
         #if KB > 176
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8352(pA), rB0
            addpd rB0, rC20
            #if KB > 177
               movddup 4120(pB), rB0
            #elif KB == 177
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8352(pA), rB1
            addpd rB1, rC21
            #if KB > 177
               movddup 4128(pB), rB1
            #elif KB == 177
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 177
               movaps 8368(pA), rA0
            #elif KB == 177
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 177
               movaps 8384(pA), rA1
            #elif KB == 177
               prefA(-128(pfA))
            #endif

            mulpd 8352(pA), rB2
            addpd rB2, rC22
            #if KB > 177
               movddup 4136(pB), rB2
            #elif KB == 177
               prefA((pfA))
            #endif

         #endif
         #if KB > 177
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8400(pA), rB0
            addpd rB0, rC20
            #if KB > 178
               movddup 4144(pB), rB0
            #elif KB == 178
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8400(pA), rB1
            addpd rB1, rC21
            #if KB > 178
               movddup 4152(pB), rB1
            #elif KB == 178
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 178
               movaps 8416(pA), rA0
            #elif KB == 178
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 178
               movaps 8432(pA), rA1
            #elif KB == 178
               prefA(-128(pfA))
            #endif

            mulpd 8400(pA), rB2
            addpd rB2, rC22
            #if KB > 178
               movddup 4160(pB), rB2
            #elif KB == 178
               prefA((pfA))
            #endif

         #endif
         #if KB > 178
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8448(pA), rB0
            addpd rB0, rC20
            #if KB > 179
               movddup 4168(pB), rB0
            #elif KB == 179
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8448(pA), rB1
            addpd rB1, rC21
            #if KB > 179
               movddup 4176(pB), rB1
            #elif KB == 179
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 179
               movaps 8464(pA), rA0
            #elif KB == 179
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 179
               movaps 8480(pA), rA1
            #elif KB == 179
               prefA(-128(pfA))
            #endif

            mulpd 8448(pA), rB2
            addpd rB2, rC22
            #if KB > 179
               movddup 4184(pB), rB2
            #elif KB == 179
               prefA((pfA))
            #endif

         #endif
         #if KB > 179
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8496(pA), rB0
            addpd rB0, rC20
            #if KB > 180
               movddup 4192(pB), rB0
            #elif KB == 180
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8496(pA), rB1
            addpd rB1, rC21
            #if KB > 180
               movddup 4200(pB), rB1
            #elif KB == 180
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 180
               movaps 8512(pA), rA0
            #elif KB == 180
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 180
               movaps 8528(pA), rA1
            #elif KB == 180
               prefA(-128(pfA))
            #endif

            mulpd 8496(pA), rB2
            addpd rB2, rC22
            #if KB > 180
               movddup 4208(pB), rB2
            #elif KB == 180
               prefA((pfA))
            #endif

         #endif
         #if KB > 180
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8544(pA), rB0
            addpd rB0, rC20
            #if KB > 181
               movddup 4216(pB), rB0
            #elif KB == 181
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8544(pA), rB1
            addpd rB1, rC21
            #if KB > 181
               movddup 4224(pB), rB1
            #elif KB == 181
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 181
               movaps 8560(pA), rA0
            #elif KB == 181
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 181
               movaps 8576(pA), rA1
            #elif KB == 181
               prefA(-128(pfA))
            #endif

            mulpd 8544(pA), rB2
            addpd rB2, rC22
            #if KB > 181
               movddup 4232(pB), rB2
            #elif KB == 181
               prefA((pfA))
            #endif

         #endif
         #if KB > 181
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8592(pA), rB0
            addpd rB0, rC20
            #if KB > 182
               movddup 4240(pB), rB0
            #elif KB == 182
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8592(pA), rB1
            addpd rB1, rC21
            #if KB > 182
               movddup 4248(pB), rB1
            #elif KB == 182
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 182
               movaps 8608(pA), rA0
            #elif KB == 182
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 182
               movaps 8624(pA), rA1
            #elif KB == 182
               prefA(-128(pfA))
            #endif

            mulpd 8592(pA), rB2
            addpd rB2, rC22
            #if KB > 182
               movddup 4256(pB), rB2
            #elif KB == 182
               prefA((pfA))
            #endif

         #endif
         #if KB > 182
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8640(pA), rB0
            addpd rB0, rC20
            #if KB > 183
               movddup 4264(pB), rB0
            #elif KB == 183
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8640(pA), rB1
            addpd rB1, rC21
            #if KB > 183
               movddup 4272(pB), rB1
            #elif KB == 183
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 183
               movaps 8656(pA), rA0
            #elif KB == 183
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 183
               movaps 8672(pA), rA1
            #elif KB == 183
               prefA(-128(pfA))
            #endif

            mulpd 8640(pA), rB2
            addpd rB2, rC22
            #if KB > 183
               movddup 4280(pB), rB2
            #elif KB == 183
               prefA((pfA))
            #endif

         #endif
         #if KB > 183
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8688(pA), rB0
            addpd rB0, rC20
            #if KB > 184
               movddup 4288(pB), rB0
            #elif KB == 184
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8688(pA), rB1
            addpd rB1, rC21
            #if KB > 184
               movddup 4296(pB), rB1
            #elif KB == 184
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 184
               movaps 8704(pA), rA0
            #elif KB == 184
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 184
               movaps 8720(pA), rA1
            #elif KB == 184
               prefA(-128(pfA))
            #endif

            mulpd 8688(pA), rB2
            addpd rB2, rC22
            #if KB > 184
               movddup 4304(pB), rB2
            #elif KB == 184
               prefA((pfA))
            #endif

         #endif
         #if KB > 184
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8736(pA), rB0
            addpd rB0, rC20
            #if KB > 185
               movddup 4312(pB), rB0
            #elif KB == 185
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8736(pA), rB1
            addpd rB1, rC21
            #if KB > 185
               movddup 4320(pB), rB1
            #elif KB == 185
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 185
               movaps 8752(pA), rA0
            #elif KB == 185
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 185
               movaps 8768(pA), rA1
            #elif KB == 185
               prefA(-128(pfA))
            #endif

            mulpd 8736(pA), rB2
            addpd rB2, rC22
            #if KB > 185
               movddup 4328(pB), rB2
            #elif KB == 185
               prefA((pfA))
            #endif

         #endif
         #if KB > 185
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8784(pA), rB0
            addpd rB0, rC20
            #if KB > 186
               movddup 4336(pB), rB0
            #elif KB == 186
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8784(pA), rB1
            addpd rB1, rC21
            #if KB > 186
               movddup 4344(pB), rB1
            #elif KB == 186
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 186
               movaps 8800(pA), rA0
            #elif KB == 186
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 186
               movaps 8816(pA), rA1
            #elif KB == 186
               prefA(-128(pfA))
            #endif

            mulpd 8784(pA), rB2
            addpd rB2, rC22
            #if KB > 186
               movddup 4352(pB), rB2
            #elif KB == 186
               prefA((pfA))
            #endif

         #endif
         #if KB > 186
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8832(pA), rB0
            addpd rB0, rC20
            #if KB > 187
               movddup 4360(pB), rB0
            #elif KB == 187
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8832(pA), rB1
            addpd rB1, rC21
            #if KB > 187
               movddup 4368(pB), rB1
            #elif KB == 187
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 187
               movaps 8848(pA), rA0
            #elif KB == 187
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 187
               movaps 8864(pA), rA1
            #elif KB == 187
               prefA(-128(pfA))
            #endif

            mulpd 8832(pA), rB2
            addpd rB2, rC22
            #if KB > 187
               movddup 4376(pB), rB2
            #elif KB == 187
               prefA((pfA))
            #endif

         #endif
         #if KB > 187
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8880(pA), rB0
            addpd rB0, rC20
            #if KB > 188
               movddup 4384(pB), rB0
            #elif KB == 188
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8880(pA), rB1
            addpd rB1, rC21
            #if KB > 188
               movddup 4392(pB), rB1
            #elif KB == 188
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 188
               movaps 8896(pA), rA0
            #elif KB == 188
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 188
               movaps 8912(pA), rA1
            #elif KB == 188
               prefA(-128(pfA))
            #endif

            mulpd 8880(pA), rB2
            addpd rB2, rC22
            #if KB > 188
               movddup 4400(pB), rB2
            #elif KB == 188
               prefA((pfA))
            #endif

         #endif
         #if KB > 188
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8928(pA), rB0
            addpd rB0, rC20
            #if KB > 189
               movddup 4408(pB), rB0
            #elif KB == 189
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8928(pA), rB1
            addpd rB1, rC21
            #if KB > 189
               movddup 4416(pB), rB1
            #elif KB == 189
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 189
               movaps 8944(pA), rA0
            #elif KB == 189
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 189
               movaps 8960(pA), rA1
            #elif KB == 189
               prefA(-128(pfA))
            #endif

            mulpd 8928(pA), rB2
            addpd rB2, rC22
            #if KB > 189
               movddup 4424(pB), rB2
            #elif KB == 189
               prefA((pfA))
            #endif

         #endif
         #if KB > 189
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 8976(pA), rB0
            addpd rB0, rC20
            #if KB > 190
               movddup 4432(pB), rB0
            #elif KB == 190
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 8976(pA), rB1
            addpd rB1, rC21
            #if KB > 190
               movddup 4440(pB), rB1
            #elif KB == 190
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 190
               movaps 8992(pA), rA0
            #elif KB == 190
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 190
               movaps 9008(pA), rA1
            #elif KB == 190
               prefA(-128(pfA))
            #endif

            mulpd 8976(pA), rB2
            addpd rB2, rC22
            #if KB > 190
               movddup 4448(pB), rB2
            #elif KB == 190
               prefA((pfA))
            #endif

         #endif
         #if KB > 190
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9024(pA), rB0
            addpd rB0, rC20
            #if KB > 191
               movddup 4456(pB), rB0
            #elif KB == 191
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9024(pA), rB1
            addpd rB1, rC21
            #if KB > 191
               movddup 4464(pB), rB1
            #elif KB == 191
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 191
               movaps 9040(pA), rA0
            #elif KB == 191
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 191
               movaps 9056(pA), rA1
            #elif KB == 191
               prefA(-128(pfA))
            #endif

            mulpd 9024(pA), rB2
            addpd rB2, rC22
            #if KB > 191
               movddup 4472(pB), rB2
            #elif KB == 191
               prefA((pfA))
            #endif

         #endif
         #if KB > 191
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9072(pA), rB0
            addpd rB0, rC20
            #if KB > 192
               movddup 4480(pB), rB0
            #elif KB == 192
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9072(pA), rB1
            addpd rB1, rC21
            #if KB > 192
               movddup 4488(pB), rB1
            #elif KB == 192
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 192
               movaps 9088(pA), rA0
            #elif KB == 192
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 192
               movaps 9104(pA), rA1
            #elif KB == 192
               prefA(-128(pfA))
            #endif

            mulpd 9072(pA), rB2
            addpd rB2, rC22
            #if KB > 192
               movddup 4496(pB), rB2
            #elif KB == 192
               prefA((pfA))
            #endif

         #endif
         #if KB > 192
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9120(pA), rB0
            addpd rB0, rC20
            #if KB > 193
               movddup 4504(pB), rB0
            #elif KB == 193
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9120(pA), rB1
            addpd rB1, rC21
            #if KB > 193
               movddup 4512(pB), rB1
            #elif KB == 193
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 193
               movaps 9136(pA), rA0
            #elif KB == 193
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 193
               movaps 9152(pA), rA1
            #elif KB == 193
               prefA(-128(pfA))
            #endif

            mulpd 9120(pA), rB2
            addpd rB2, rC22
            #if KB > 193
               movddup 4520(pB), rB2
            #elif KB == 193
               prefA((pfA))
            #endif

         #endif
         #if KB > 193
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9168(pA), rB0
            addpd rB0, rC20
            #if KB > 194
               movddup 4528(pB), rB0
            #elif KB == 194
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9168(pA), rB1
            addpd rB1, rC21
            #if KB > 194
               movddup 4536(pB), rB1
            #elif KB == 194
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 194
               movaps 9184(pA), rA0
            #elif KB == 194
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 194
               movaps 9200(pA), rA1
            #elif KB == 194
               prefA(-128(pfA))
            #endif

            mulpd 9168(pA), rB2
            addpd rB2, rC22
            #if KB > 194
               movddup 4544(pB), rB2
            #elif KB == 194
               prefA((pfA))
            #endif

         #endif
         #if KB > 194
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9216(pA), rB0
            addpd rB0, rC20
            #if KB > 195
               movddup 4552(pB), rB0
            #elif KB == 195
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9216(pA), rB1
            addpd rB1, rC21
            #if KB > 195
               movddup 4560(pB), rB1
            #elif KB == 195
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 195
               movaps 9232(pA), rA0
            #elif KB == 195
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 195
               movaps 9248(pA), rA1
            #elif KB == 195
               prefA(-128(pfA))
            #endif

            mulpd 9216(pA), rB2
            addpd rB2, rC22
            #if KB > 195
               movddup 4568(pB), rB2
            #elif KB == 195
               prefA((pfA))
            #endif

         #endif
         #if KB > 195
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9264(pA), rB0
            addpd rB0, rC20
            #if KB > 196
               movddup 4576(pB), rB0
            #elif KB == 196
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9264(pA), rB1
            addpd rB1, rC21
            #if KB > 196
               movddup 4584(pB), rB1
            #elif KB == 196
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 196
               movaps 9280(pA), rA0
            #elif KB == 196
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 196
               movaps 9296(pA), rA1
            #elif KB == 196
               prefA(-128(pfA))
            #endif

            mulpd 9264(pA), rB2
            addpd rB2, rC22
            #if KB > 196
               movddup 4592(pB), rB2
            #elif KB == 196
               prefA((pfA))
            #endif

         #endif
         #if KB > 196
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9312(pA), rB0
            addpd rB0, rC20
            #if KB > 197
               movddup 4600(pB), rB0
            #elif KB == 197
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9312(pA), rB1
            addpd rB1, rC21
            #if KB > 197
               movddup 4608(pB), rB1
            #elif KB == 197
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 197
               movaps 9328(pA), rA0
            #elif KB == 197
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 197
               movaps 9344(pA), rA1
            #elif KB == 197
               prefA(-128(pfA))
            #endif

            mulpd 9312(pA), rB2
            addpd rB2, rC22
            #if KB > 197
               movddup 4616(pB), rB2
            #elif KB == 197
               prefA((pfA))
            #endif

         #endif
         #if KB > 197
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9360(pA), rB0
            addpd rB0, rC20
            #if KB > 198
               movddup 4624(pB), rB0
            #elif KB == 198
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9360(pA), rB1
            addpd rB1, rC21
            #if KB > 198
               movddup 4632(pB), rB1
            #elif KB == 198
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 198
               movaps 9376(pA), rA0
            #elif KB == 198
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 198
               movaps 9392(pA), rA1
            #elif KB == 198
               prefA(-128(pfA))
            #endif

            mulpd 9360(pA), rB2
            addpd rB2, rC22
            #if KB > 198
               movddup 4640(pB), rB2
            #elif KB == 198
               prefA((pfA))
            #endif

         #endif
         #if KB > 198
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9408(pA), rB0
            addpd rB0, rC20
            #if KB > 199
               movddup 4648(pB), rB0
            #elif KB == 199
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9408(pA), rB1
            addpd rB1, rC21
            #if KB > 199
               movddup 4656(pB), rB1
            #elif KB == 199
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 199
               movaps 9424(pA), rA0
            #elif KB == 199
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 199
               movaps 9440(pA), rA1
            #elif KB == 199
               prefA(-128(pfA))
            #endif

            mulpd 9408(pA), rB2
            addpd rB2, rC22
            #if KB > 199
               movddup 4664(pB), rB2
            #elif KB == 199
               prefA((pfA))
            #endif

         #endif
         #if KB > 199
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9456(pA), rB0
            addpd rB0, rC20
            #if KB > 200
               movddup 4672(pB), rB0
            #elif KB == 200
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9456(pA), rB1
            addpd rB1, rC21
            #if KB > 200
               movddup 4680(pB), rB1
            #elif KB == 200
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 200
               movaps 9472(pA), rA0
            #elif KB == 200
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 200
               movaps 9488(pA), rA1
            #elif KB == 200
               prefA(-128(pfA))
            #endif

            mulpd 9456(pA), rB2
            addpd rB2, rC22
            #if KB > 200
               movddup 4688(pB), rB2
            #elif KB == 200
               prefA((pfA))
            #endif

         #endif
         #if KB > 200
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9504(pA), rB0
            addpd rB0, rC20
            #if KB > 201
               movddup 4696(pB), rB0
            #elif KB == 201
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9504(pA), rB1
            addpd rB1, rC21
            #if KB > 201
               movddup 4704(pB), rB1
            #elif KB == 201
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 201
               movaps 9520(pA), rA0
            #elif KB == 201
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 201
               movaps 9536(pA), rA1
            #elif KB == 201
               prefA(-128(pfA))
            #endif

            mulpd 9504(pA), rB2
            addpd rB2, rC22
            #if KB > 201
               movddup 4712(pB), rB2
            #elif KB == 201
               prefA((pfA))
            #endif

         #endif
         #if KB > 201
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9552(pA), rB0
            addpd rB0, rC20
            #if KB > 202
               movddup 4720(pB), rB0
            #elif KB == 202
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9552(pA), rB1
            addpd rB1, rC21
            #if KB > 202
               movddup 4728(pB), rB1
            #elif KB == 202
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 202
               movaps 9568(pA), rA0
            #elif KB == 202
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 202
               movaps 9584(pA), rA1
            #elif KB == 202
               prefA(-128(pfA))
            #endif

            mulpd 9552(pA), rB2
            addpd rB2, rC22
            #if KB > 202
               movddup 4736(pB), rB2
            #elif KB == 202
               prefA((pfA))
            #endif

         #endif
         #if KB > 202
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9600(pA), rB0
            addpd rB0, rC20
            #if KB > 203
               movddup 4744(pB), rB0
            #elif KB == 203
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9600(pA), rB1
            addpd rB1, rC21
            #if KB > 203
               movddup 4752(pB), rB1
            #elif KB == 203
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 203
               movaps 9616(pA), rA0
            #elif KB == 203
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 203
               movaps 9632(pA), rA1
            #elif KB == 203
               prefA(-128(pfA))
            #endif

            mulpd 9600(pA), rB2
            addpd rB2, rC22
            #if KB > 203
               movddup 4760(pB), rB2
            #elif KB == 203
               prefA((pfA))
            #endif

         #endif
         #if KB > 203
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9648(pA), rB0
            addpd rB0, rC20
            #if KB > 204
               movddup 4768(pB), rB0
            #elif KB == 204
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9648(pA), rB1
            addpd rB1, rC21
            #if KB > 204
               movddup 4776(pB), rB1
            #elif KB == 204
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 204
               movaps 9664(pA), rA0
            #elif KB == 204
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 204
               movaps 9680(pA), rA1
            #elif KB == 204
               prefA(-128(pfA))
            #endif

            mulpd 9648(pA), rB2
            addpd rB2, rC22
            #if KB > 204
               movddup 4784(pB), rB2
            #elif KB == 204
               prefA((pfA))
            #endif

         #endif
         #if KB > 204
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9696(pA), rB0
            addpd rB0, rC20
            #if KB > 205
               movddup 4792(pB), rB0
            #elif KB == 205
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9696(pA), rB1
            addpd rB1, rC21
            #if KB > 205
               movddup 4800(pB), rB1
            #elif KB == 205
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 205
               movaps 9712(pA), rA0
            #elif KB == 205
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 205
               movaps 9728(pA), rA1
            #elif KB == 205
               prefA(-128(pfA))
            #endif

            mulpd 9696(pA), rB2
            addpd rB2, rC22
            #if KB > 205
               movddup 4808(pB), rB2
            #elif KB == 205
               prefA((pfA))
            #endif

         #endif
         #if KB > 205
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9744(pA), rB0
            addpd rB0, rC20
            #if KB > 206
               movddup 4816(pB), rB0
            #elif KB == 206
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9744(pA), rB1
            addpd rB1, rC21
            #if KB > 206
               movddup 4824(pB), rB1
            #elif KB == 206
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 206
               movaps 9760(pA), rA0
            #elif KB == 206
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 206
               movaps 9776(pA), rA1
            #elif KB == 206
               prefA(-128(pfA))
            #endif

            mulpd 9744(pA), rB2
            addpd rB2, rC22
            #if KB > 206
               movddup 4832(pB), rB2
            #elif KB == 206
               prefA((pfA))
            #endif

         #endif
         #if KB > 206
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9792(pA), rB0
            addpd rB0, rC20
            #if KB > 207
               movddup 4840(pB), rB0
            #elif KB == 207
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9792(pA), rB1
            addpd rB1, rC21
            #if KB > 207
               movddup 4848(pB), rB1
            #elif KB == 207
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 207
               movaps 9808(pA), rA0
            #elif KB == 207
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 207
               movaps 9824(pA), rA1
            #elif KB == 207
               prefA(-128(pfA))
            #endif

            mulpd 9792(pA), rB2
            addpd rB2, rC22
            #if KB > 207
               movddup 4856(pB), rB2
            #elif KB == 207
               prefA((pfA))
            #endif

         #endif
         #if KB > 207
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9840(pA), rB0
            addpd rB0, rC20
            #if KB > 208
               movddup 4864(pB), rB0
            #elif KB == 208
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9840(pA), rB1
            addpd rB1, rC21
            #if KB > 208
               movddup 4872(pB), rB1
            #elif KB == 208
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 208
               movaps 9856(pA), rA0
            #elif KB == 208
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 208
               movaps 9872(pA), rA1
            #elif KB == 208
               prefA(-128(pfA))
            #endif

            mulpd 9840(pA), rB2
            addpd rB2, rC22
            #if KB > 208
               movddup 4880(pB), rB2
            #elif KB == 208
               prefA((pfA))
            #endif

         #endif
         #if KB > 208
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9888(pA), rB0
            addpd rB0, rC20
            #if KB > 209
               movddup 4888(pB), rB0
            #elif KB == 209
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9888(pA), rB1
            addpd rB1, rC21
            #if KB > 209
               movddup 4896(pB), rB1
            #elif KB == 209
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 209
               movaps 9904(pA), rA0
            #elif KB == 209
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 209
               movaps 9920(pA), rA1
            #elif KB == 209
               prefA(-128(pfA))
            #endif

            mulpd 9888(pA), rB2
            addpd rB2, rC22
            #if KB > 209
               movddup 4904(pB), rB2
            #elif KB == 209
               prefA((pfA))
            #endif

         #endif
         #if KB > 209
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9936(pA), rB0
            addpd rB0, rC20
            #if KB > 210
               movddup 4912(pB), rB0
            #elif KB == 210
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9936(pA), rB1
            addpd rB1, rC21
            #if KB > 210
               movddup 4920(pB), rB1
            #elif KB == 210
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 210
               movaps 9952(pA), rA0
            #elif KB == 210
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 210
               movaps 9968(pA), rA1
            #elif KB == 210
               prefA(-128(pfA))
            #endif

            mulpd 9936(pA), rB2
            addpd rB2, rC22
            #if KB > 210
               movddup 4928(pB), rB2
            #elif KB == 210
               prefA((pfA))
            #endif

         #endif
         #if KB > 210
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 9984(pA), rB0
            addpd rB0, rC20
            #if KB > 211
               movddup 4936(pB), rB0
            #elif KB == 211
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 9984(pA), rB1
            addpd rB1, rC21
            #if KB > 211
               movddup 4944(pB), rB1
            #elif KB == 211
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 211
               movaps 10000(pA), rA0
            #elif KB == 211
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 211
               movaps 10016(pA), rA1
            #elif KB == 211
               prefA(-128(pfA))
            #endif

            mulpd 9984(pA), rB2
            addpd rB2, rC22
            #if KB > 211
               movddup 4952(pB), rB2
            #elif KB == 211
               prefA((pfA))
            #endif

         #endif
         #if KB > 211
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10032(pA), rB0
            addpd rB0, rC20
            #if KB > 212
               movddup 4960(pB), rB0
            #elif KB == 212
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10032(pA), rB1
            addpd rB1, rC21
            #if KB > 212
               movddup 4968(pB), rB1
            #elif KB == 212
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 212
               movaps 10048(pA), rA0
            #elif KB == 212
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 212
               movaps 10064(pA), rA1
            #elif KB == 212
               prefA(-128(pfA))
            #endif

            mulpd 10032(pA), rB2
            addpd rB2, rC22
            #if KB > 212
               movddup 4976(pB), rB2
            #elif KB == 212
               prefA((pfA))
            #endif

         #endif
         #if KB > 212
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10080(pA), rB0
            addpd rB0, rC20
            #if KB > 213
               movddup 4984(pB), rB0
            #elif KB == 213
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10080(pA), rB1
            addpd rB1, rC21
            #if KB > 213
               movddup 4992(pB), rB1
            #elif KB == 213
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 213
               movaps 10096(pA), rA0
            #elif KB == 213
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 213
               movaps 10112(pA), rA1
            #elif KB == 213
               prefA(-128(pfA))
            #endif

            mulpd 10080(pA), rB2
            addpd rB2, rC22
            #if KB > 213
               movddup 5000(pB), rB2
            #elif KB == 213
               prefA((pfA))
            #endif

         #endif
         #if KB > 213
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10128(pA), rB0
            addpd rB0, rC20
            #if KB > 214
               movddup 5008(pB), rB0
            #elif KB == 214
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10128(pA), rB1
            addpd rB1, rC21
            #if KB > 214
               movddup 5016(pB), rB1
            #elif KB == 214
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 214
               movaps 10144(pA), rA0
            #elif KB == 214
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 214
               movaps 10160(pA), rA1
            #elif KB == 214
               prefA(-128(pfA))
            #endif

            mulpd 10128(pA), rB2
            addpd rB2, rC22
            #if KB > 214
               movddup 5024(pB), rB2
            #elif KB == 214
               prefA((pfA))
            #endif

         #endif
         #if KB > 214
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10176(pA), rB0
            addpd rB0, rC20
            #if KB > 215
               movddup 5032(pB), rB0
            #elif KB == 215
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10176(pA), rB1
            addpd rB1, rC21
            #if KB > 215
               movddup 5040(pB), rB1
            #elif KB == 215
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 215
               movaps 10192(pA), rA0
            #elif KB == 215
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 215
               movaps 10208(pA), rA1
            #elif KB == 215
               prefA(-128(pfA))
            #endif

            mulpd 10176(pA), rB2
            addpd rB2, rC22
            #if KB > 215
               movddup 5048(pB), rB2
            #elif KB == 215
               prefA((pfA))
            #endif

         #endif
         #if KB > 215
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10224(pA), rB0
            addpd rB0, rC20
            #if KB > 216
               movddup 5056(pB), rB0
            #elif KB == 216
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10224(pA), rB1
            addpd rB1, rC21
            #if KB > 216
               movddup 5064(pB), rB1
            #elif KB == 216
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 216
               movaps 10240(pA), rA0
            #elif KB == 216
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 216
               movaps 10256(pA), rA1
            #elif KB == 216
               prefA(-128(pfA))
            #endif

            mulpd 10224(pA), rB2
            addpd rB2, rC22
            #if KB > 216
               movddup 5072(pB), rB2
            #elif KB == 216
               prefA((pfA))
            #endif

         #endif
         #if KB > 216
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10272(pA), rB0
            addpd rB0, rC20
            #if KB > 217
               movddup 5080(pB), rB0
            #elif KB == 217
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10272(pA), rB1
            addpd rB1, rC21
            #if KB > 217
               movddup 5088(pB), rB1
            #elif KB == 217
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 217
               movaps 10288(pA), rA0
            #elif KB == 217
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 217
               movaps 10304(pA), rA1
            #elif KB == 217
               prefA(-128(pfA))
            #endif

            mulpd 10272(pA), rB2
            addpd rB2, rC22
            #if KB > 217
               movddup 5096(pB), rB2
            #elif KB == 217
               prefA((pfA))
            #endif

         #endif
         #if KB > 217
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10320(pA), rB0
            addpd rB0, rC20
            #if KB > 218
               movddup 5104(pB), rB0
            #elif KB == 218
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10320(pA), rB1
            addpd rB1, rC21
            #if KB > 218
               movddup 5112(pB), rB1
            #elif KB == 218
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 218
               movaps 10336(pA), rA0
            #elif KB == 218
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 218
               movaps 10352(pA), rA1
            #elif KB == 218
               prefA(-128(pfA))
            #endif

            mulpd 10320(pA), rB2
            addpd rB2, rC22
            #if KB > 218
               movddup 5120(pB), rB2
            #elif KB == 218
               prefA((pfA))
            #endif

         #endif
         #if KB > 218
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10368(pA), rB0
            addpd rB0, rC20
            #if KB > 219
               movddup 5128(pB), rB0
            #elif KB == 219
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10368(pA), rB1
            addpd rB1, rC21
            #if KB > 219
               movddup 5136(pB), rB1
            #elif KB == 219
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 219
               movaps 10384(pA), rA0
            #elif KB == 219
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 219
               movaps 10400(pA), rA1
            #elif KB == 219
               prefA(-128(pfA))
            #endif

            mulpd 10368(pA), rB2
            addpd rB2, rC22
            #if KB > 219
               movddup 5144(pB), rB2
            #elif KB == 219
               prefA((pfA))
            #endif

         #endif
         #if KB > 219
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10416(pA), rB0
            addpd rB0, rC20
            #if KB > 220
               movddup 5152(pB), rB0
            #elif KB == 220
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10416(pA), rB1
            addpd rB1, rC21
            #if KB > 220
               movddup 5160(pB), rB1
            #elif KB == 220
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 220
               movaps 10432(pA), rA0
            #elif KB == 220
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 220
               movaps 10448(pA), rA1
            #elif KB == 220
               prefA(-128(pfA))
            #endif

            mulpd 10416(pA), rB2
            addpd rB2, rC22
            #if KB > 220
               movddup 5168(pB), rB2
            #elif KB == 220
               prefA((pfA))
            #endif

         #endif
         #if KB > 220
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10464(pA), rB0
            addpd rB0, rC20
            #if KB > 221
               movddup 5176(pB), rB0
            #elif KB == 221
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10464(pA), rB1
            addpd rB1, rC21
            #if KB > 221
               movddup 5184(pB), rB1
            #elif KB == 221
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 221
               movaps 10480(pA), rA0
            #elif KB == 221
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 221
               movaps 10496(pA), rA1
            #elif KB == 221
               prefA(-128(pfA))
            #endif

            mulpd 10464(pA), rB2
            addpd rB2, rC22
            #if KB > 221
               movddup 5192(pB), rB2
            #elif KB == 221
               prefA((pfA))
            #endif

         #endif
         #if KB > 221
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10512(pA), rB0
            addpd rB0, rC20
            #if KB > 222
               movddup 5200(pB), rB0
            #elif KB == 222
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10512(pA), rB1
            addpd rB1, rC21
            #if KB > 222
               movddup 5208(pB), rB1
            #elif KB == 222
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 222
               movaps 10528(pA), rA0
            #elif KB == 222
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 222
               movaps 10544(pA), rA1
            #elif KB == 222
               prefA(-128(pfA))
            #endif

            mulpd 10512(pA), rB2
            addpd rB2, rC22
            #if KB > 222
               movddup 5216(pB), rB2
            #elif KB == 222
               prefA((pfA))
            #endif

         #endif
         #if KB > 222
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10560(pA), rB0
            addpd rB0, rC20
            #if KB > 223
               movddup 5224(pB), rB0
            #elif KB == 223
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10560(pA), rB1
            addpd rB1, rC21
            #if KB > 223
               movddup 5232(pB), rB1
            #elif KB == 223
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 223
               movaps 10576(pA), rA0
            #elif KB == 223
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 223
               movaps 10592(pA), rA1
            #elif KB == 223
               prefA(-128(pfA))
            #endif

            mulpd 10560(pA), rB2
            addpd rB2, rC22
            #if KB > 223
               movddup 5240(pB), rB2
            #elif KB == 223
               prefA((pfA))
            #endif

         #endif
         #if KB > 223
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10608(pA), rB0
            addpd rB0, rC20
            #if KB > 224
               movddup 5248(pB), rB0
            #elif KB == 224
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10608(pA), rB1
            addpd rB1, rC21
            #if KB > 224
               movddup 5256(pB), rB1
            #elif KB == 224
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 224
               movaps 10624(pA), rA0
            #elif KB == 224
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 224
               movaps 10640(pA), rA1
            #elif KB == 224
               prefA(-128(pfA))
            #endif

            mulpd 10608(pA), rB2
            addpd rB2, rC22
            #if KB > 224
               movddup 5264(pB), rB2
            #elif KB == 224
               prefA((pfA))
            #endif

         #endif
         #if KB > 224
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10656(pA), rB0
            addpd rB0, rC20
            #if KB > 225
               movddup 5272(pB), rB0
            #elif KB == 225
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10656(pA), rB1
            addpd rB1, rC21
            #if KB > 225
               movddup 5280(pB), rB1
            #elif KB == 225
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 225
               movaps 10672(pA), rA0
            #elif KB == 225
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 225
               movaps 10688(pA), rA1
            #elif KB == 225
               prefA(-128(pfA))
            #endif

            mulpd 10656(pA), rB2
            addpd rB2, rC22
            #if KB > 225
               movddup 5288(pB), rB2
            #elif KB == 225
               prefA((pfA))
            #endif

         #endif
         #if KB > 225
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10704(pA), rB0
            addpd rB0, rC20
            #if KB > 226
               movddup 5296(pB), rB0
            #elif KB == 226
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10704(pA), rB1
            addpd rB1, rC21
            #if KB > 226
               movddup 5304(pB), rB1
            #elif KB == 226
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 226
               movaps 10720(pA), rA0
            #elif KB == 226
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 226
               movaps 10736(pA), rA1
            #elif KB == 226
               prefA(-128(pfA))
            #endif

            mulpd 10704(pA), rB2
            addpd rB2, rC22
            #if KB > 226
               movddup 5312(pB), rB2
            #elif KB == 226
               prefA((pfA))
            #endif

         #endif
         #if KB > 226
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10752(pA), rB0
            addpd rB0, rC20
            #if KB > 227
               movddup 5320(pB), rB0
            #elif KB == 227
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10752(pA), rB1
            addpd rB1, rC21
            #if KB > 227
               movddup 5328(pB), rB1
            #elif KB == 227
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 227
               movaps 10768(pA), rA0
            #elif KB == 227
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 227
               movaps 10784(pA), rA1
            #elif KB == 227
               prefA(-128(pfA))
            #endif

            mulpd 10752(pA), rB2
            addpd rB2, rC22
            #if KB > 227
               movddup 5336(pB), rB2
            #elif KB == 227
               prefA((pfA))
            #endif

         #endif
         #if KB > 227
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10800(pA), rB0
            addpd rB0, rC20
            #if KB > 228
               movddup 5344(pB), rB0
            #elif KB == 228
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10800(pA), rB1
            addpd rB1, rC21
            #if KB > 228
               movddup 5352(pB), rB1
            #elif KB == 228
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 228
               movaps 10816(pA), rA0
            #elif KB == 228
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 228
               movaps 10832(pA), rA1
            #elif KB == 228
               prefA(-128(pfA))
            #endif

            mulpd 10800(pA), rB2
            addpd rB2, rC22
            #if KB > 228
               movddup 5360(pB), rB2
            #elif KB == 228
               prefA((pfA))
            #endif

         #endif
         #if KB > 228
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10848(pA), rB0
            addpd rB0, rC20
            #if KB > 229
               movddup 5368(pB), rB0
            #elif KB == 229
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10848(pA), rB1
            addpd rB1, rC21
            #if KB > 229
               movddup 5376(pB), rB1
            #elif KB == 229
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 229
               movaps 10864(pA), rA0
            #elif KB == 229
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 229
               movaps 10880(pA), rA1
            #elif KB == 229
               prefA(-128(pfA))
            #endif

            mulpd 10848(pA), rB2
            addpd rB2, rC22
            #if KB > 229
               movddup 5384(pB), rB2
            #elif KB == 229
               prefA((pfA))
            #endif

         #endif
         #if KB > 229
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10896(pA), rB0
            addpd rB0, rC20
            #if KB > 230
               movddup 5392(pB), rB0
            #elif KB == 230
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10896(pA), rB1
            addpd rB1, rC21
            #if KB > 230
               movddup 5400(pB), rB1
            #elif KB == 230
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 230
               movaps 10912(pA), rA0
            #elif KB == 230
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 230
               movaps 10928(pA), rA1
            #elif KB == 230
               prefA(-128(pfA))
            #endif

            mulpd 10896(pA), rB2
            addpd rB2, rC22
            #if KB > 230
               movddup 5408(pB), rB2
            #elif KB == 230
               prefA((pfA))
            #endif

         #endif
         #if KB > 230
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10944(pA), rB0
            addpd rB0, rC20
            #if KB > 231
               movddup 5416(pB), rB0
            #elif KB == 231
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10944(pA), rB1
            addpd rB1, rC21
            #if KB > 231
               movddup 5424(pB), rB1
            #elif KB == 231
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 231
               movaps 10960(pA), rA0
            #elif KB == 231
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 231
               movaps 10976(pA), rA1
            #elif KB == 231
               prefA(-128(pfA))
            #endif

            mulpd 10944(pA), rB2
            addpd rB2, rC22
            #if KB > 231
               movddup 5432(pB), rB2
            #elif KB == 231
               prefA((pfA))
            #endif

         #endif
         #if KB > 231
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 10992(pA), rB0
            addpd rB0, rC20
            #if KB > 232
               movddup 5440(pB), rB0
            #elif KB == 232
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 10992(pA), rB1
            addpd rB1, rC21
            #if KB > 232
               movddup 5448(pB), rB1
            #elif KB == 232
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 232
               movaps 11008(pA), rA0
            #elif KB == 232
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 232
               movaps 11024(pA), rA1
            #elif KB == 232
               prefA(-128(pfA))
            #endif

            mulpd 10992(pA), rB2
            addpd rB2, rC22
            #if KB > 232
               movddup 5456(pB), rB2
            #elif KB == 232
               prefA((pfA))
            #endif

         #endif
         #if KB > 232
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11040(pA), rB0
            addpd rB0, rC20
            #if KB > 233
               movddup 5464(pB), rB0
            #elif KB == 233
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11040(pA), rB1
            addpd rB1, rC21
            #if KB > 233
               movddup 5472(pB), rB1
            #elif KB == 233
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 233
               movaps 11056(pA), rA0
            #elif KB == 233
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 233
               movaps 11072(pA), rA1
            #elif KB == 233
               prefA(-128(pfA))
            #endif

            mulpd 11040(pA), rB2
            addpd rB2, rC22
            #if KB > 233
               movddup 5480(pB), rB2
            #elif KB == 233
               prefA((pfA))
            #endif

         #endif
         #if KB > 233
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11088(pA), rB0
            addpd rB0, rC20
            #if KB > 234
               movddup 5488(pB), rB0
            #elif KB == 234
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11088(pA), rB1
            addpd rB1, rC21
            #if KB > 234
               movddup 5496(pB), rB1
            #elif KB == 234
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 234
               movaps 11104(pA), rA0
            #elif KB == 234
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 234
               movaps 11120(pA), rA1
            #elif KB == 234
               prefA(-128(pfA))
            #endif

            mulpd 11088(pA), rB2
            addpd rB2, rC22
            #if KB > 234
               movddup 5504(pB), rB2
            #elif KB == 234
               prefA((pfA))
            #endif

         #endif
         #if KB > 234
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11136(pA), rB0
            addpd rB0, rC20
            #if KB > 235
               movddup 5512(pB), rB0
            #elif KB == 235
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11136(pA), rB1
            addpd rB1, rC21
            #if KB > 235
               movddup 5520(pB), rB1
            #elif KB == 235
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 235
               movaps 11152(pA), rA0
            #elif KB == 235
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 235
               movaps 11168(pA), rA1
            #elif KB == 235
               prefA(-128(pfA))
            #endif

            mulpd 11136(pA), rB2
            addpd rB2, rC22
            #if KB > 235
               movddup 5528(pB), rB2
            #elif KB == 235
               prefA((pfA))
            #endif

         #endif
         #if KB > 235
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11184(pA), rB0
            addpd rB0, rC20
            #if KB > 236
               movddup 5536(pB), rB0
            #elif KB == 236
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11184(pA), rB1
            addpd rB1, rC21
            #if KB > 236
               movddup 5544(pB), rB1
            #elif KB == 236
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 236
               movaps 11200(pA), rA0
            #elif KB == 236
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 236
               movaps 11216(pA), rA1
            #elif KB == 236
               prefA(-128(pfA))
            #endif

            mulpd 11184(pA), rB2
            addpd rB2, rC22
            #if KB > 236
               movddup 5552(pB), rB2
            #elif KB == 236
               prefA((pfA))
            #endif

         #endif
         #if KB > 236
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11232(pA), rB0
            addpd rB0, rC20
            #if KB > 237
               movddup 5560(pB), rB0
            #elif KB == 237
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11232(pA), rB1
            addpd rB1, rC21
            #if KB > 237
               movddup 5568(pB), rB1
            #elif KB == 237
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 237
               movaps 11248(pA), rA0
            #elif KB == 237
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 237
               movaps 11264(pA), rA1
            #elif KB == 237
               prefA(-128(pfA))
            #endif

            mulpd 11232(pA), rB2
            addpd rB2, rC22
            #if KB > 237
               movddup 5576(pB), rB2
            #elif KB == 237
               prefA((pfA))
            #endif

         #endif
         #if KB > 237
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11280(pA), rB0
            addpd rB0, rC20
            #if KB > 238
               movddup 5584(pB), rB0
            #elif KB == 238
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11280(pA), rB1
            addpd rB1, rC21
            #if KB > 238
               movddup 5592(pB), rB1
            #elif KB == 238
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 238
               movaps 11296(pA), rA0
            #elif KB == 238
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 238
               movaps 11312(pA), rA1
            #elif KB == 238
               prefA(-128(pfA))
            #endif

            mulpd 11280(pA), rB2
            addpd rB2, rC22
            #if KB > 238
               movddup 5600(pB), rB2
            #elif KB == 238
               prefA((pfA))
            #endif

         #endif
         #if KB > 238
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11328(pA), rB0
            addpd rB0, rC20
            #if KB > 239
               movddup 5608(pB), rB0
            #elif KB == 239
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11328(pA), rB1
            addpd rB1, rC21
            #if KB > 239
               movddup 5616(pB), rB1
            #elif KB == 239
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 239
               movaps 11344(pA), rA0
            #elif KB == 239
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 239
               movaps 11360(pA), rA1
            #elif KB == 239
               prefA(-128(pfA))
            #endif

            mulpd 11328(pA), rB2
            addpd rB2, rC22
            #if KB > 239
               movddup 5624(pB), rB2
            #elif KB == 239
               prefA((pfA))
            #endif

         #endif
         #if KB > 239
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11376(pA), rB0
            addpd rB0, rC20
            #if KB > 240
               movddup 5632(pB), rB0
            #elif KB == 240
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11376(pA), rB1
            addpd rB1, rC21
            #if KB > 240
               movddup 5640(pB), rB1
            #elif KB == 240
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 240
               movaps 11392(pA), rA0
            #elif KB == 240
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 240
               movaps 11408(pA), rA1
            #elif KB == 240
               prefA(-128(pfA))
            #endif

            mulpd 11376(pA), rB2
            addpd rB2, rC22
            #if KB > 240
               movddup 5648(pB), rB2
            #elif KB == 240
               prefA((pfA))
            #endif

         #endif
         #if KB > 240
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11424(pA), rB0
            addpd rB0, rC20
            #if KB > 241
               movddup 5656(pB), rB0
            #elif KB == 241
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11424(pA), rB1
            addpd rB1, rC21
            #if KB > 241
               movddup 5664(pB), rB1
            #elif KB == 241
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 241
               movaps 11440(pA), rA0
            #elif KB == 241
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 241
               movaps 11456(pA), rA1
            #elif KB == 241
               prefA(-128(pfA))
            #endif

            mulpd 11424(pA), rB2
            addpd rB2, rC22
            #if KB > 241
               movddup 5672(pB), rB2
            #elif KB == 241
               prefA((pfA))
            #endif

         #endif
         #if KB > 241
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11472(pA), rB0
            addpd rB0, rC20
            #if KB > 242
               movddup 5680(pB), rB0
            #elif KB == 242
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11472(pA), rB1
            addpd rB1, rC21
            #if KB > 242
               movddup 5688(pB), rB1
            #elif KB == 242
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 242
               movaps 11488(pA), rA0
            #elif KB == 242
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 242
               movaps 11504(pA), rA1
            #elif KB == 242
               prefA(-128(pfA))
            #endif

            mulpd 11472(pA), rB2
            addpd rB2, rC22
            #if KB > 242
               movddup 5696(pB), rB2
            #elif KB == 242
               prefA((pfA))
            #endif

         #endif
         #if KB > 242
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11520(pA), rB0
            addpd rB0, rC20
            #if KB > 243
               movddup 5704(pB), rB0
            #elif KB == 243
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11520(pA), rB1
            addpd rB1, rC21
            #if KB > 243
               movddup 5712(pB), rB1
            #elif KB == 243
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 243
               movaps 11536(pA), rA0
            #elif KB == 243
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 243
               movaps 11552(pA), rA1
            #elif KB == 243
               prefA(-128(pfA))
            #endif

            mulpd 11520(pA), rB2
            addpd rB2, rC22
            #if KB > 243
               movddup 5720(pB), rB2
            #elif KB == 243
               prefA((pfA))
            #endif

         #endif
         #if KB > 243
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11568(pA), rB0
            addpd rB0, rC20
            #if KB > 244
               movddup 5728(pB), rB0
            #elif KB == 244
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11568(pA), rB1
            addpd rB1, rC21
            #if KB > 244
               movddup 5736(pB), rB1
            #elif KB == 244
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 244
               movaps 11584(pA), rA0
            #elif KB == 244
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 244
               movaps 11600(pA), rA1
            #elif KB == 244
               prefA(-128(pfA))
            #endif

            mulpd 11568(pA), rB2
            addpd rB2, rC22
            #if KB > 244
               movddup 5744(pB), rB2
            #elif KB == 244
               prefA((pfA))
            #endif

         #endif
         #if KB > 244
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11616(pA), rB0
            addpd rB0, rC20
            #if KB > 245
               movddup 5752(pB), rB0
            #elif KB == 245
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11616(pA), rB1
            addpd rB1, rC21
            #if KB > 245
               movddup 5760(pB), rB1
            #elif KB == 245
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 245
               movaps 11632(pA), rA0
            #elif KB == 245
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 245
               movaps 11648(pA), rA1
            #elif KB == 245
               prefA(-128(pfA))
            #endif

            mulpd 11616(pA), rB2
            addpd rB2, rC22
            #if KB > 245
               movddup 5768(pB), rB2
            #elif KB == 245
               prefA((pfA))
            #endif

         #endif
         #if KB > 245
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11664(pA), rB0
            addpd rB0, rC20
            #if KB > 246
               movddup 5776(pB), rB0
            #elif KB == 246
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11664(pA), rB1
            addpd rB1, rC21
            #if KB > 246
               movddup 5784(pB), rB1
            #elif KB == 246
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 246
               movaps 11680(pA), rA0
            #elif KB == 246
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 246
               movaps 11696(pA), rA1
            #elif KB == 246
               prefA(-128(pfA))
            #endif

            mulpd 11664(pA), rB2
            addpd rB2, rC22
            #if KB > 246
               movddup 5792(pB), rB2
            #elif KB == 246
               prefA((pfA))
            #endif

         #endif
         #if KB > 246
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11712(pA), rB0
            addpd rB0, rC20
            #if KB > 247
               movddup 5800(pB), rB0
            #elif KB == 247
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11712(pA), rB1
            addpd rB1, rC21
            #if KB > 247
               movddup 5808(pB), rB1
            #elif KB == 247
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 247
               movaps 11728(pA), rA0
            #elif KB == 247
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 247
               movaps 11744(pA), rA1
            #elif KB == 247
               prefA(-128(pfA))
            #endif

            mulpd 11712(pA), rB2
            addpd rB2, rC22
            #if KB > 247
               movddup 5816(pB), rB2
            #elif KB == 247
               prefA((pfA))
            #endif

         #endif
         #if KB > 247
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11760(pA), rB0
            addpd rB0, rC20
            #if KB > 248
               movddup 5824(pB), rB0
            #elif KB == 248
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11760(pA), rB1
            addpd rB1, rC21
            #if KB > 248
               movddup 5832(pB), rB1
            #elif KB == 248
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 248
               movaps 11776(pA), rA0
            #elif KB == 248
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 248
               movaps 11792(pA), rA1
            #elif KB == 248
               prefA(-128(pfA))
            #endif

            mulpd 11760(pA), rB2
            addpd rB2, rC22
            #if KB > 248
               movddup 5840(pB), rB2
            #elif KB == 248
               prefA((pfA))
            #endif

         #endif
         #if KB > 248
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11808(pA), rB0
            addpd rB0, rC20
            #if KB > 249
               movddup 5848(pB), rB0
            #elif KB == 249
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11808(pA), rB1
            addpd rB1, rC21
            #if KB > 249
               movddup 5856(pB), rB1
            #elif KB == 249
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 249
               movaps 11824(pA), rA0
            #elif KB == 249
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 249
               movaps 11840(pA), rA1
            #elif KB == 249
               prefA(-128(pfA))
            #endif

            mulpd 11808(pA), rB2
            addpd rB2, rC22
            #if KB > 249
               movddup 5864(pB), rB2
            #elif KB == 249
               prefA((pfA))
            #endif

         #endif
         #if KB > 249
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11856(pA), rB0
            addpd rB0, rC20
            #if KB > 250
               movddup 5872(pB), rB0
            #elif KB == 250
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11856(pA), rB1
            addpd rB1, rC21
            #if KB > 250
               movddup 5880(pB), rB1
            #elif KB == 250
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 250
               movaps 11872(pA), rA0
            #elif KB == 250
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 250
               movaps 11888(pA), rA1
            #elif KB == 250
               prefA(-128(pfA))
            #endif

            mulpd 11856(pA), rB2
            addpd rB2, rC22
            #if KB > 250
               movddup 5888(pB), rB2
            #elif KB == 250
               prefA((pfA))
            #endif

         #endif
         #if KB > 250
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11904(pA), rB0
            addpd rB0, rC20
            #if KB > 251
               movddup 5896(pB), rB0
            #elif KB == 251
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11904(pA), rB1
            addpd rB1, rC21
            #if KB > 251
               movddup 5904(pB), rB1
            #elif KB == 251
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 251
               movaps 11920(pA), rA0
            #elif KB == 251
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 251
               movaps 11936(pA), rA1
            #elif KB == 251
               prefA(-128(pfA))
            #endif

            mulpd 11904(pA), rB2
            addpd rB2, rC22
            #if KB > 251
               movddup 5912(pB), rB2
            #elif KB == 251
               prefA((pfA))
            #endif

         #endif
         #if KB > 251
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 11952(pA), rB0
            addpd rB0, rC20
            #if KB > 252
               movddup 5920(pB), rB0
            #elif KB == 252
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 11952(pA), rB1
            addpd rB1, rC21
            #if KB > 252
               movddup 5928(pB), rB1
            #elif KB == 252
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 252
               movaps 11968(pA), rA0
            #elif KB == 252
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 252
               movaps 11984(pA), rA1
            #elif KB == 252
               prefA(-128(pfA))
            #endif

            mulpd 11952(pA), rB2
            addpd rB2, rC22
            #if KB > 252
               movddup 5936(pB), rB2
            #elif KB == 252
               prefA((pfA))
            #endif

         #endif
         #if KB > 252
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 12000(pA), rB0
            addpd rB0, rC20
            #if KB > 253
               movddup 5944(pB), rB0
            #elif KB == 253
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 12000(pA), rB1
            addpd rB1, rC21
            #if KB > 253
               movddup 5952(pB), rB1
            #elif KB == 253
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 253
               movaps 12016(pA), rA0
            #elif KB == 253
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 253
               movaps 12032(pA), rA1
            #elif KB == 253
               prefA(-128(pfA))
            #endif

            mulpd 12000(pA), rB2
            addpd rB2, rC22
            #if KB > 253
               movddup 5960(pB), rB2
            #elif KB == 253
               prefA((pfA))
            #endif

         #endif
         #if KB > 253
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 12048(pA), rB0
            addpd rB0, rC20
            #if KB > 254
               movddup 5968(pB), rB0
            #elif KB == 254
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 12048(pA), rB1
            addpd rB1, rC21
            #if KB > 254
               movddup 5976(pB), rB1
            #elif KB == 254
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 254
               movaps 12064(pA), rA0
            #elif KB == 254
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 254
               movaps 12080(pA), rA1
            #elif KB == 254
               prefA(-128(pfA))
            #endif

            mulpd 12048(pA), rB2
            addpd rB2, rC22
            #if KB > 254
               movddup 5984(pB), rB2
            #elif KB == 254
               prefA((pfA))
            #endif

         #endif
         #if KB > 254
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 12096(pA), rB0
            addpd rB0, rC20
            #if KB > 255
               movddup 5992(pB), rB0
            #elif KB == 255
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 12096(pA), rB1
            addpd rB1, rC21
            #if KB > 255
               movddup 6000(pB), rB1
            #elif KB == 255
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 255
               movaps 12112(pA), rA0
            #elif KB == 255
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 255
               movaps 12128(pA), rA1
            #elif KB == 255
               prefA(-128(pfA))
            #endif

            mulpd 12096(pA), rB2
            addpd rB2, rC22
            #if KB > 255
               movddup 6008(pB), rB2
            #elif KB == 255
               prefA((pfA))
            #endif

         #endif
         #if KB > 255
            movapd rA0, rm0
            mulpd rB0, rm0
            addpd rm0, rC00

            movapd rA1, rm0
            mulpd rB0, rm0
            addpd rm0, rC10

            mulpd 12144(pA), rB0
            addpd rB0, rC20
            #if KB > 256
               movddup 6016(pB), rB0
            #elif KB == 256
               prefB(-128(pfB))
            #endif

            movapd rA0, rm0
            mulpd rB1, rm0
            addpd rm0, rC01

            movapd rA1, rm0
            mulpd rB1, rm0
            addpd rm0, rC11

            mulpd 12144(pA), rB1
            addpd rB1, rC21
            #if KB > 256
               movddup 6024(pB), rB1
            #elif KB == 256
               prefB(-64(pfB))
            #endif

            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 256
               movaps 12160(pA), rA0
            #elif KB == 256
               prefB((pfB))
            #endif

            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 256
               movaps 12176(pA), rA1
            #elif KB == 256
               prefA(-128(pfA))
            #endif

            mulpd 12144(pA), rB2
            addpd rB2, rC22
            #if KB > 256
               movddup 6032(pB), rB2
            #elif KB == 256
               prefA((pfA))
            #endif

         #endif

         add incPF, pfA
         add incPF, pfB
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
