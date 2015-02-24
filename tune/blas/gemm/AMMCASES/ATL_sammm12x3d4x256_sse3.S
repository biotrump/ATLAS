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
         movapd rC01, rC11
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
               prefetcht0 -48+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps -48(pB), rB2
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
               movaps -32(pB), rB0
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
               movaps -16(pB), rB1
            #endif

         #endif
         #if KB > 2
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 0(pA), rA2
               prefetcht0 0(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 0(pB), rB2
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
               movaps 16(pB), rB0
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
               movaps 32(pB), rB1
            #endif

         #endif
         #if KB > 3
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 48(pA), rA2
               prefetcht0 48+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 48(pB), rB2
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
               movaps 64(pB), rB0
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
               movaps 80(pB), rB1
            #endif

         #endif
         #if KB > 4
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 96(pA), rA2
               prefetcht0 96(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 96(pB), rB2
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
               movaps 112(pB), rB0
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
               movaps 128(pB), rB1
            #endif

         #endif
         #if KB > 5
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 144(pA), rA2
               prefetcht0 144+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 144(pB), rB2
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
               movaps 160(pB), rB0
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
               movaps 176(pB), rB1
            #endif

         #endif
         #if KB > 6
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 192(pA), rA2
               prefetcht0 192(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 192(pB), rB2
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
               movaps 208(pB), rB0
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
               movaps 224(pB), rB1
            #endif

         #endif
         #if KB > 7
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 240(pA), rA2
               prefetcht0 240+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 240(pB), rB2
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
               movaps 256(pB), rB0
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
               movaps 272(pB), rB1
            #endif

         #endif
         #if KB > 8
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 288(pA), rA2
               prefetcht0 288(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 288(pB), rB2
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
               movaps 304(pB), rB0
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
               movaps 320(pB), rB1
            #endif

         #endif
         #if KB > 9
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 336(pA), rA2
               prefetcht0 336+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 336(pB), rB2
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
               movaps 352(pB), rB0
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
               movaps 368(pB), rB1
            #endif

         #endif
         #if KB > 10
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 384(pA), rA2
               prefetcht0 384(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 384(pB), rB2
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
               movaps 400(pB), rB0
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
               movaps 416(pB), rB1
            #endif

         #endif
         #if KB > 11
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 432(pA), rA2
               prefetcht0 432+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 432(pB), rB2
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
               movaps 448(pB), rB0
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
               movaps 464(pB), rB1
            #endif

         #endif
         #if KB > 12
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 480(pA), rA2
               prefetcht0 480(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 480(pB), rB2
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
               movaps 496(pB), rB0
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
               movaps 512(pB), rB1
            #endif

         #endif
         #if KB > 13
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 528(pA), rA2
               prefetcht0 528+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 528(pB), rB2
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
               movaps 544(pB), rB0
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
               movaps 560(pB), rB1
            #endif

         #endif
         #if KB > 14
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 576(pA), rA2
               prefetcht0 576(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 576(pB), rB2
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
               movaps 592(pB), rB0
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
               movaps 608(pB), rB1
            #endif

         #endif
         #if KB > 15
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 624(pA), rA2
               prefetcht0 624+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 624(pB), rB2
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
               movaps 640(pB), rB0
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
               movaps 656(pB), rB1
            #endif

         #endif
         #if KB > 16
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 672(pA), rA2
               prefetcht0 672(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 672(pB), rB2
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
               movaps 688(pB), rB0
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
               movaps 704(pB), rB1
            #endif

         #endif
         #if KB > 17
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 720(pA), rA2
               prefetcht0 720+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 720(pB), rB2
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
               movaps 736(pB), rB0
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
               movaps 752(pB), rB1
            #endif

         #endif
         #if KB > 18
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 768(pA), rA2
               prefetcht0 768(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 768(pB), rB2
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
               movaps 784(pB), rB0
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
               movaps 800(pB), rB1
            #endif

         #endif
         #if KB > 19
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 816(pA), rA2
               prefetcht0 816+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 816(pB), rB2
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
               movaps 832(pB), rB0
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
               movaps 848(pB), rB1
            #endif

         #endif
         #if KB > 20
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 864(pA), rA2
               prefetcht0 864(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 864(pB), rB2
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
               movaps 880(pB), rB0
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
               movaps 896(pB), rB1
            #endif

         #endif
         #if KB > 21
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 912(pA), rA2
               prefetcht0 912+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 912(pB), rB2
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
               movaps 928(pB), rB0
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
               movaps 944(pB), rB1
            #endif

         #endif
         #if KB > 22
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 960(pA), rA2
               prefetcht0 960(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 960(pB), rB2
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
               movaps 976(pB), rB0
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
               movaps 992(pB), rB1
            #endif

         #endif
         #if KB > 23
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1008(pA), rA2
               prefetcht0 1008+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1008(pB), rB2
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
               movaps 1024(pB), rB0
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
               movaps 1040(pB), rB1
            #endif

         #endif
         #if KB > 24
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1056(pA), rA2
               prefetcht0 1056(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1056(pB), rB2
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
               movaps 1072(pB), rB0
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
               movaps 1088(pB), rB1
            #endif

         #endif
         #if KB > 25
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1104(pA), rA2
               prefetcht0 1104+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1104(pB), rB2
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
               movaps 1120(pB), rB0
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
               movaps 1136(pB), rB1
            #endif

         #endif
         #if KB > 26
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1152(pA), rA2
               prefetcht0 1152(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1152(pB), rB2
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
               movaps 1168(pB), rB0
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
               movaps 1184(pB), rB1
            #endif

         #endif
         #if KB > 27
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1200(pA), rA2
               prefetcht0 1200+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1200(pB), rB2
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
               movaps 1216(pB), rB0
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
               movaps 1232(pB), rB1
            #endif

         #endif
         #if KB > 28
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1248(pA), rA2
               prefetcht0 1248(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1248(pB), rB2
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
               movaps 1264(pB), rB0
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
               movaps 1280(pB), rB1
            #endif

         #endif
         #if KB > 29
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1296(pA), rA2
               prefetcht0 1296+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1296(pB), rB2
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
               movaps 1312(pB), rB0
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
               movaps 1328(pB), rB1
            #endif

         #endif
         #if KB > 30
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1344(pA), rA2
               prefetcht0 1344(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1344(pB), rB2
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
               movaps 1360(pB), rB0
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
               movaps 1376(pB), rB1
            #endif

         #endif
         #if KB > 31
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1392(pA), rA2
               prefetcht0 1392+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1392(pB), rB2
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
               movaps 1408(pB), rB0
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
               movaps 1424(pB), rB1
            #endif

         #endif
         #if KB > 32
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1440(pA), rA2
               prefetcht0 1440(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1440(pB), rB2
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
               movaps 1456(pB), rB0
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
               movaps 1472(pB), rB1
            #endif

         #endif
         #if KB > 33
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1488(pA), rA2
               prefetcht0 1488+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1488(pB), rB2
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
               movaps 1504(pB), rB0
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
               movaps 1520(pB), rB1
            #endif

         #endif
         #if KB > 34
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1536(pA), rA2
               prefetcht0 1536(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1536(pB), rB2
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
               movaps 1552(pB), rB0
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
               movaps 1568(pB), rB1
            #endif

         #endif
         #if KB > 35
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1584(pA), rA2
               prefetcht0 1584+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1584(pB), rB2
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
               movaps 1600(pB), rB0
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
               movaps 1616(pB), rB1
            #endif

         #endif
         #if KB > 36
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1632(pA), rA2
               prefetcht0 1632(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1632(pB), rB2
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
               movaps 1648(pB), rB0
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
               movaps 1664(pB), rB1
            #endif

         #endif
         #if KB > 37
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1680(pA), rA2
               prefetcht0 1680+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1680(pB), rB2
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
               movaps 1696(pB), rB0
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
               movaps 1712(pB), rB1
            #endif

         #endif
         #if KB > 38
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1728(pA), rA2
               prefetcht0 1728(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1728(pB), rB2
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
               movaps 1744(pB), rB0
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
               movaps 1760(pB), rB1
            #endif

         #endif
         #if KB > 39
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1776(pA), rA2
               prefetcht0 1776+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1776(pB), rB2
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
               movaps 1792(pB), rB0
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
               movaps 1808(pB), rB1
            #endif

         #endif
         #if KB > 40
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1824(pA), rA2
               prefetcht0 1824(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1824(pB), rB2
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
               movaps 1840(pB), rB0
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
               movaps 1856(pB), rB1
            #endif

         #endif
         #if KB > 41
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1872(pA), rA2
               prefetcht0 1872+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1872(pB), rB2
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
               movaps 1888(pB), rB0
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
               movaps 1904(pB), rB1
            #endif

         #endif
         #if KB > 42
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1920(pA), rA2
               prefetcht0 1920(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1920(pB), rB2
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
               movaps 1936(pB), rB0
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
               movaps 1952(pB), rB1
            #endif

         #endif
         #if KB > 43
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 1968(pA), rA2
               prefetcht0 1968+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 1968(pB), rB2
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
               movaps 1984(pB), rB0
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
               movaps 2000(pB), rB1
            #endif

         #endif
         #if KB > 44
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2016(pA), rA2
               prefetcht0 2016(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2016(pB), rB2
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
               movaps 2032(pB), rB0
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
               movaps 2048(pB), rB1
            #endif

         #endif
         #if KB > 45
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2064(pA), rA2
               prefetcht0 2064+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2064(pB), rB2
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
               movaps 2080(pB), rB0
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
               movaps 2096(pB), rB1
            #endif

         #endif
         #if KB > 46
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2112(pA), rA2
               prefetcht0 2112(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2112(pB), rB2
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
               movaps 2128(pB), rB0
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
               movaps 2144(pB), rB1
            #endif

         #endif
         #if KB > 47
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2160(pA), rA2
               prefetcht0 2160+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2160(pB), rB2
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
               movaps 2176(pB), rB0
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
               movaps 2192(pB), rB1
            #endif

         #endif
         #if KB > 48
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2208(pA), rA2
               prefetcht0 2208(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2208(pB), rB2
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
               movaps 2224(pB), rB0
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
               movaps 2240(pB), rB1
            #endif

         #endif
         #if KB > 49
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2256(pA), rA2
               prefetcht0 2256+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2256(pB), rB2
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
               movaps 2272(pB), rB0
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
               movaps 2288(pB), rB1
            #endif

         #endif
         #if KB > 50
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2304(pA), rA2
               prefetcht0 2304(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2304(pB), rB2
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
               movaps 2320(pB), rB0
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
               movaps 2336(pB), rB1
            #endif

         #endif
         #if KB > 51
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2352(pA), rA2
               prefetcht0 2352+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2352(pB), rB2
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
               movaps 2368(pB), rB0
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
               movaps 2384(pB), rB1
            #endif

         #endif
         #if KB > 52
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2400(pA), rA2
               prefetcht0 2400(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2400(pB), rB2
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
               movaps 2416(pB), rB0
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
               movaps 2432(pB), rB1
            #endif

         #endif
         #if KB > 53
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2448(pA), rA2
               prefetcht0 2448+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2448(pB), rB2
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
               movaps 2464(pB), rB0
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
               movaps 2480(pB), rB1
            #endif

         #endif
         #if KB > 54
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2496(pA), rA2
               prefetcht0 2496(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2496(pB), rB2
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
               movaps 2512(pB), rB0
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
               movaps 2528(pB), rB1
            #endif

         #endif
         #if KB > 55
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2544(pA), rA2
               prefetcht0 2544+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2544(pB), rB2
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
               movaps 2560(pB), rB0
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
               movaps 2576(pB), rB1
            #endif

         #endif
         #if KB > 56
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2592(pA), rA2
               prefetcht0 2592(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2592(pB), rB2
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
               movaps 2608(pB), rB0
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
               movaps 2624(pB), rB1
            #endif

         #endif
         #if KB > 57
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2640(pA), rA2
               prefetcht0 2640+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2640(pB), rB2
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
               movaps 2656(pB), rB0
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
               movaps 2672(pB), rB1
            #endif

         #endif
         #if KB > 58
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2688(pA), rA2
               prefetcht0 2688(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2688(pB), rB2
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
               movaps 2704(pB), rB0
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
               movaps 2720(pB), rB1
            #endif

         #endif
         #if KB > 59
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2736(pA), rA2
               prefetcht0 2736+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2736(pB), rB2
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
               movaps 2752(pB), rB0
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
               movaps 2768(pB), rB1
            #endif

         #endif
         #if KB > 60
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2784(pA), rA2
               prefetcht0 2784(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2784(pB), rB2
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
               movaps 2800(pB), rB0
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
               movaps 2816(pB), rB1
            #endif

         #endif
         #if KB > 61
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2832(pA), rA2
               prefetcht0 2832+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2832(pB), rB2
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
               movaps 2848(pB), rB0
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
               movaps 2864(pB), rB1
            #endif

         #endif
         #if KB > 62
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2880(pA), rA2
               prefetcht0 2880(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2880(pB), rB2
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
               movaps 2896(pB), rB0
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
               movaps 2912(pB), rB1
            #endif

         #endif
         #if KB > 63
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2928(pA), rA2
               prefetcht0 2928+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2928(pB), rB2
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
               movaps 2944(pB), rB0
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
               movaps 2960(pB), rB1
            #endif

         #endif
         #if KB > 64
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 2976(pA), rA2
               prefetcht0 2976(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 2976(pB), rB2
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
               movaps 2992(pB), rB0
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
               movaps 3008(pB), rB1
            #endif

         #endif
         #if KB > 65
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3024(pA), rA2
               prefetcht0 3024+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3024(pB), rB2
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
               movaps 3040(pB), rB0
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
               movaps 3056(pB), rB1
            #endif

         #endif
         #if KB > 66
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3072(pA), rA2
               prefetcht0 3072(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3072(pB), rB2
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
               movaps 3088(pB), rB0
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
               movaps 3104(pB), rB1
            #endif

         #endif
         #if KB > 67
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3120(pA), rA2
               prefetcht0 3120+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3120(pB), rB2
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
               movaps 3136(pB), rB0
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
               movaps 3152(pB), rB1
            #endif

         #endif
         #if KB > 68
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3168(pA), rA2
               prefetcht0 3168(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3168(pB), rB2
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
               movaps 3184(pB), rB0
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
               movaps 3200(pB), rB1
            #endif

         #endif
         #if KB > 69
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3216(pA), rA2
               prefetcht0 3216+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3216(pB), rB2
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
               movaps 3232(pB), rB0
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
               movaps 3248(pB), rB1
            #endif

         #endif
         #if KB > 70
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3264(pA), rA2
               prefetcht0 3264(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3264(pB), rB2
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
               movaps 3280(pB), rB0
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
               movaps 3296(pB), rB1
            #endif

         #endif
         #if KB > 71
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3312(pA), rA2
               prefetcht0 3312+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3312(pB), rB2
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
               movaps 3328(pB), rB0
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
               movaps 3344(pB), rB1
            #endif

         #endif
         #if KB > 72
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3360(pA), rA2
               prefetcht0 3360(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3360(pB), rB2
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
               movaps 3376(pB), rB0
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
               movaps 3392(pB), rB1
            #endif

         #endif
         #if KB > 73
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3408(pA), rA2
               prefetcht0 3408+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3408(pB), rB2
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
               movaps 3424(pB), rB0
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
               movaps 3440(pB), rB1
            #endif

         #endif
         #if KB > 74
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3456(pA), rA2
               prefetcht0 3456(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3456(pB), rB2
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
               movaps 3472(pB), rB0
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
               movaps 3488(pB), rB1
            #endif

         #endif
         #if KB > 75
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3504(pA), rA2
               prefetcht0 3504+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3504(pB), rB2
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
               movaps 3520(pB), rB0
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
               movaps 3536(pB), rB1
            #endif

         #endif
         #if KB > 76
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3552(pA), rA2
               prefetcht0 3552(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3552(pB), rB2
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
               movaps 3568(pB), rB0
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
               movaps 3584(pB), rB1
            #endif

         #endif
         #if KB > 77
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3600(pA), rA2
               prefetcht0 3600+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3600(pB), rB2
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
               movaps 3616(pB), rB0
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
               movaps 3632(pB), rB1
            #endif

         #endif
         #if KB > 78
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3648(pA), rA2
               prefetcht0 3648(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3648(pB), rB2
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
               movaps 3664(pB), rB0
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
               movaps 3680(pB), rB1
            #endif

         #endif
         #if KB > 79
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3696(pA), rA2
               prefetcht0 3696+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3696(pB), rB2
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
               movaps 3712(pB), rB0
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
               movaps 3728(pB), rB1
            #endif

         #endif
         #if KB > 80
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3744(pA), rA2
               prefetcht0 3744(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3744(pB), rB2
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
               movaps 3760(pB), rB0
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
               movaps 3776(pB), rB1
            #endif

         #endif
         #if KB > 81
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3792(pA), rA2
               prefetcht0 3792+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3792(pB), rB2
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
               movaps 3808(pB), rB0
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
               movaps 3824(pB), rB1
            #endif

         #endif
         #if KB > 82
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3840(pA), rA2
               prefetcht0 3840(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3840(pB), rB2
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
               movaps 3856(pB), rB0
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
               movaps 3872(pB), rB1
            #endif

         #endif
         #if KB > 83
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3888(pA), rA2
               prefetcht0 3888+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3888(pB), rB2
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
               movaps 3904(pB), rB0
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
               movaps 3920(pB), rB1
            #endif

         #endif
         #if KB > 84
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3936(pA), rA2
               prefetcht0 3936(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3936(pB), rB2
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
               movaps 3952(pB), rB0
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
               movaps 3968(pB), rB1
            #endif

         #endif
         #if KB > 85
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 3984(pA), rA2
               prefetcht0 3984+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 3984(pB), rB2
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
               movaps 4000(pB), rB0
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
               movaps 4016(pB), rB1
            #endif

         #endif
         #if KB > 86
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4032(pA), rA2
               prefetcht0 4032(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4032(pB), rB2
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
               movaps 4048(pB), rB0
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
               movaps 4064(pB), rB1
            #endif

         #endif
         #if KB > 87
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4080(pA), rA2
               prefetcht0 4080+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4080(pB), rB2
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
               movaps 4096(pB), rB0
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
               movaps 4112(pB), rB1
            #endif

         #endif
         #if KB > 88
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4128(pA), rA2
               prefetcht0 4128(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4128(pB), rB2
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
               movaps 4144(pB), rB0
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
               movaps 4160(pB), rB1
            #endif

         #endif
         #if KB > 89
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4176(pA), rA2
               prefetcht0 4176+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4176(pB), rB2
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
               movaps 4192(pB), rB0
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
               movaps 4208(pB), rB1
            #endif

         #endif
         #if KB > 90
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4224(pA), rA2
               prefetcht0 4224(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4224(pB), rB2
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
               movaps 4240(pB), rB0
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
               movaps 4256(pB), rB1
            #endif

         #endif
         #if KB > 91
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4272(pA), rA2
               prefetcht0 4272+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4272(pB), rB2
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
               movaps 4288(pB), rB0
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
               movaps 4304(pB), rB1
            #endif

         #endif
         #if KB > 92
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4320(pA), rA2
               prefetcht0 4320(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4320(pB), rB2
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
               movaps 4336(pB), rB0
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
               movaps 4352(pB), rB1
            #endif

         #endif
         #if KB > 93
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4368(pA), rA2
               prefetcht0 4368+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4368(pB), rB2
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
               movaps 4384(pB), rB0
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
               movaps 4400(pB), rB1
            #endif

         #endif
         #if KB > 94
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4416(pA), rA2
               prefetcht0 4416(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4416(pB), rB2
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
               movaps 4432(pB), rB0
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
               movaps 4448(pB), rB1
            #endif

         #endif
         #if KB > 95
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4464(pA), rA2
               prefetcht0 4464+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4464(pB), rB2
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
               movaps 4480(pB), rB0
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
               movaps 4496(pB), rB1
            #endif

         #endif
         #if KB > 96
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4512(pA), rA2
               prefetcht0 4512(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4512(pB), rB2
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
               movaps 4528(pB), rB0
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
               movaps 4544(pB), rB1
            #endif

         #endif
         #if KB > 97
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4560(pA), rA2
               prefetcht0 4560+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4560(pB), rB2
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
               movaps 4576(pB), rB0
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
               movaps 4592(pB), rB1
            #endif

         #endif
         #if KB > 98
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4608(pA), rA2
               prefetcht0 4608(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4608(pB), rB2
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
               movaps 4624(pB), rB0
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
               movaps 4640(pB), rB1
            #endif

         #endif
         #if KB > 99
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4656(pA), rA2
               prefetcht0 4656+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4656(pB), rB2
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
               movaps 4672(pB), rB0
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
               movaps 4688(pB), rB1
            #endif

         #endif
         #if KB > 100
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4704(pA), rA2
               prefetcht0 4704(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4704(pB), rB2
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
               movaps 4720(pB), rB0
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
               movaps 4736(pB), rB1
            #endif

         #endif
         #if KB > 101
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4752(pA), rA2
               prefetcht0 4752+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4752(pB), rB2
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
               movaps 4768(pB), rB0
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
               movaps 4784(pB), rB1
            #endif

         #endif
         #if KB > 102
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4800(pA), rA2
               prefetcht0 4800(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4800(pB), rB2
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
               movaps 4816(pB), rB0
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
               movaps 4832(pB), rB1
            #endif

         #endif
         #if KB > 103
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4848(pA), rA2
               prefetcht0 4848+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4848(pB), rB2
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
               movaps 4864(pB), rB0
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
               movaps 4880(pB), rB1
            #endif

         #endif
         #if KB > 104
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4896(pA), rA2
               prefetcht0 4896(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4896(pB), rB2
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
               movaps 4912(pB), rB0
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
               movaps 4928(pB), rB1
            #endif

         #endif
         #if KB > 105
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4944(pA), rA2
               prefetcht0 4944+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4944(pB), rB2
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
               movaps 4960(pB), rB0
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
               movaps 4976(pB), rB1
            #endif

         #endif
         #if KB > 106
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 4992(pA), rA2
               prefetcht0 4992(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 4992(pB), rB2
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
               movaps 5008(pB), rB0
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
               movaps 5024(pB), rB1
            #endif

         #endif
         #if KB > 107
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5040(pA), rA2
               prefetcht0 5040+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5040(pB), rB2
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
               movaps 5056(pB), rB0
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
               movaps 5072(pB), rB1
            #endif

         #endif
         #if KB > 108
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5088(pA), rA2
               prefetcht0 5088(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5088(pB), rB2
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
               movaps 5104(pB), rB0
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
               movaps 5120(pB), rB1
            #endif

         #endif
         #if KB > 109
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5136(pA), rA2
               prefetcht0 5136+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5136(pB), rB2
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
               movaps 5152(pB), rB0
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
               movaps 5168(pB), rB1
            #endif

         #endif
         #if KB > 110
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5184(pA), rA2
               prefetcht0 5184(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5184(pB), rB2
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
               movaps 5200(pB), rB0
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
               movaps 5216(pB), rB1
            #endif

         #endif
         #if KB > 111
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5232(pA), rA2
               prefetcht0 5232+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5232(pB), rB2
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
               movaps 5248(pB), rB0
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
               movaps 5264(pB), rB1
            #endif

         #endif
         #if KB > 112
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5280(pA), rA2
               prefetcht0 5280(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5280(pB), rB2
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
               movaps 5296(pB), rB0
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
               movaps 5312(pB), rB1
            #endif

         #endif
         #if KB > 113
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5328(pA), rA2
               prefetcht0 5328+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5328(pB), rB2
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
               movaps 5344(pB), rB0
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
               movaps 5360(pB), rB1
            #endif

         #endif
         #if KB > 114
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5376(pA), rA2
               prefetcht0 5376(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5376(pB), rB2
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
               movaps 5392(pB), rB0
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
               movaps 5408(pB), rB1
            #endif

         #endif
         #if KB > 115
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5424(pA), rA2
               prefetcht0 5424+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5424(pB), rB2
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
               movaps 5440(pB), rB0
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
               movaps 5456(pB), rB1
            #endif

         #endif
         #if KB > 116
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5472(pA), rA2
               prefetcht0 5472(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5472(pB), rB2
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
               movaps 5488(pB), rB0
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
               movaps 5504(pB), rB1
            #endif

         #endif
         #if KB > 117
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5520(pA), rA2
               prefetcht0 5520+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5520(pB), rB2
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
               movaps 5536(pB), rB0
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
               movaps 5552(pB), rB1
            #endif

         #endif
         #if KB > 118
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5568(pA), rA2
               prefetcht0 5568(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5568(pB), rB2
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
               movaps 5584(pB), rB0
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
               movaps 5600(pB), rB1
            #endif

         #endif
         #if KB > 119
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5616(pA), rA2
               prefetcht0 5616+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5616(pB), rB2
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
               movaps 5632(pB), rB0
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
               movaps 5648(pB), rB1
            #endif

         #endif
         #if KB > 120
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5664(pA), rA2
               prefetcht0 5664(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5664(pB), rB2
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
               movaps 5680(pB), rB0
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
               movaps 5696(pB), rB1
            #endif

         #endif
         #if KB > 121
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5712(pA), rA2
               prefetcht0 5712+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5712(pB), rB2
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
               movaps 5728(pB), rB0
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
               movaps 5744(pB), rB1
            #endif

         #endif
         #if KB > 122
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5760(pA), rA2
               prefetcht0 5760(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5760(pB), rB2
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
               movaps 5776(pB), rB0
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
               movaps 5792(pB), rB1
            #endif

         #endif
         #if KB > 123
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5808(pA), rA2
               prefetcht0 5808+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5808(pB), rB2
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
               movaps 5824(pB), rB0
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
               movaps 5840(pB), rB1
            #endif

         #endif
         #if KB > 124
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5856(pA), rA2
               prefetcht0 5856(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5856(pB), rB2
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
               movaps 5872(pB), rB0
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
               movaps 5888(pB), rB1
            #endif

         #endif
         #if KB > 125
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5904(pA), rA2
               prefetcht0 5904+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5904(pB), rB2
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
               movaps 5920(pB), rB0
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
               movaps 5936(pB), rB1
            #endif

         #endif
         #if KB > 126
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 5952(pA), rA2
               prefetcht0 5952(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 5952(pB), rB2
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
               movaps 5968(pB), rB0
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
               movaps 5984(pB), rB1
            #endif

         #endif
         #if KB > 127
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6000(pA), rA2
               prefetcht0 6000+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6000(pB), rB2
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
               movaps 6016(pB), rB0
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
               movaps 6032(pB), rB1
            #endif

         #endif
         #if KB > 128
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6048(pA), rA2
               prefetcht0 6048(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6048(pB), rB2
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
               movaps 6064(pB), rB0
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
               movaps 6080(pB), rB1
            #endif

         #endif
         #if KB > 129
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6096(pA), rA2
               prefetcht0 6096+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6096(pB), rB2
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
               movaps 6112(pB), rB0
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
               movaps 6128(pB), rB1
            #endif

         #endif
         #if KB > 130
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6144(pA), rA2
               prefetcht0 6144(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6144(pB), rB2
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
               movaps 6160(pB), rB0
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
               movaps 6176(pB), rB1
            #endif

         #endif
         #if KB > 131
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6192(pA), rA2
               prefetcht0 6192+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6192(pB), rB2
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
               movaps 6208(pB), rB0
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
               movaps 6224(pB), rB1
            #endif

         #endif
         #if KB > 132
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6240(pA), rA2
               prefetcht0 6240(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6240(pB), rB2
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
               movaps 6256(pB), rB0
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
               movaps 6272(pB), rB1
            #endif

         #endif
         #if KB > 133
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6288(pA), rA2
               prefetcht0 6288+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6288(pB), rB2
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
               movaps 6304(pB), rB0
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
               movaps 6320(pB), rB1
            #endif

         #endif
         #if KB > 134
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6336(pA), rA2
               prefetcht0 6336(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6336(pB), rB2
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
               movaps 6352(pB), rB0
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
               movaps 6368(pB), rB1
            #endif

         #endif
         #if KB > 135
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6384(pA), rA2
               prefetcht0 6384+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6384(pB), rB2
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
               movaps 6400(pB), rB0
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
               movaps 6416(pB), rB1
            #endif

         #endif
         #if KB > 136
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6432(pA), rA2
               prefetcht0 6432(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6432(pB), rB2
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
               movaps 6448(pB), rB0
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
               movaps 6464(pB), rB1
            #endif

         #endif
         #if KB > 137
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6480(pA), rA2
               prefetcht0 6480+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6480(pB), rB2
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
               movaps 6496(pB), rB0
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
               movaps 6512(pB), rB1
            #endif

         #endif
         #if KB > 138
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6528(pA), rA2
               prefetcht0 6528(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6528(pB), rB2
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
               movaps 6544(pB), rB0
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
               movaps 6560(pB), rB1
            #endif

         #endif
         #if KB > 139
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6576(pA), rA2
               prefetcht0 6576+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6576(pB), rB2
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
               movaps 6592(pB), rB0
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
               movaps 6608(pB), rB1
            #endif

         #endif
         #if KB > 140
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6624(pA), rA2
               prefetcht0 6624(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6624(pB), rB2
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
               movaps 6640(pB), rB0
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
               movaps 6656(pB), rB1
            #endif

         #endif
         #if KB > 141
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6672(pA), rA2
               prefetcht0 6672+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6672(pB), rB2
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
               movaps 6688(pB), rB0
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
               movaps 6704(pB), rB1
            #endif

         #endif
         #if KB > 142
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6720(pA), rA2
               prefetcht0 6720(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6720(pB), rB2
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
               movaps 6736(pB), rB0
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
               movaps 6752(pB), rB1
            #endif

         #endif
         #if KB > 143
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6768(pA), rA2
               prefetcht0 6768+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6768(pB), rB2
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
               movaps 6784(pB), rB0
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
               movaps 6800(pB), rB1
            #endif

         #endif
         #if KB > 144
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6816(pA), rA2
               prefetcht0 6816(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6816(pB), rB2
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
               movaps 6832(pB), rB0
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
               movaps 6848(pB), rB1
            #endif

         #endif
         #if KB > 145
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6864(pA), rA2
               prefetcht0 6864+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6864(pB), rB2
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
               movaps 6880(pB), rB0
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
               movaps 6896(pB), rB1
            #endif

         #endif
         #if KB > 146
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6912(pA), rA2
               prefetcht0 6912(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6912(pB), rB2
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
               movaps 6928(pB), rB0
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
               movaps 6944(pB), rB1
            #endif

         #endif
         #if KB > 147
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 6960(pA), rA2
               prefetcht0 6960+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 6960(pB), rB2
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
               movaps 6976(pB), rB0
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
               movaps 6992(pB), rB1
            #endif

         #endif
         #if KB > 148
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7008(pA), rA2
               prefetcht0 7008(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7008(pB), rB2
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
               movaps 7024(pB), rB0
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
               movaps 7040(pB), rB1
            #endif

         #endif
         #if KB > 149
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7056(pA), rA2
               prefetcht0 7056+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7056(pB), rB2
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
               movaps 7072(pB), rB0
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
               movaps 7088(pB), rB1
            #endif

         #endif
         #if KB > 150
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7104(pA), rA2
               prefetcht0 7104(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7104(pB), rB2
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
               movaps 7120(pB), rB0
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
               movaps 7136(pB), rB1
            #endif

         #endif
         #if KB > 151
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7152(pA), rA2
               prefetcht0 7152+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7152(pB), rB2
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
               movaps 7168(pB), rB0
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
               movaps 7184(pB), rB1
            #endif

         #endif
         #if KB > 152
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7200(pA), rA2
               prefetcht0 7200(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7200(pB), rB2
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
               movaps 7216(pB), rB0
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
               movaps 7232(pB), rB1
            #endif

         #endif
         #if KB > 153
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7248(pA), rA2
               prefetcht0 7248+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7248(pB), rB2
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
               movaps 7264(pB), rB0
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
               movaps 7280(pB), rB1
            #endif

         #endif
         #if KB > 154
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7296(pA), rA2
               prefetcht0 7296(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7296(pB), rB2
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
               movaps 7312(pB), rB0
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
               movaps 7328(pB), rB1
            #endif

         #endif
         #if KB > 155
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7344(pA), rA2
               prefetcht0 7344+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7344(pB), rB2
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
               movaps 7360(pB), rB0
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
               movaps 7376(pB), rB1
            #endif

         #endif
         #if KB > 156
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7392(pA), rA2
               prefetcht0 7392(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7392(pB), rB2
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
               movaps 7408(pB), rB0
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
               movaps 7424(pB), rB1
            #endif

         #endif
         #if KB > 157
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7440(pA), rA2
               prefetcht0 7440+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7440(pB), rB2
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
               movaps 7456(pB), rB0
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
               movaps 7472(pB), rB1
            #endif

         #endif
         #if KB > 158
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7488(pA), rA2
               prefetcht0 7488(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7488(pB), rB2
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
               movaps 7504(pB), rB0
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
               movaps 7520(pB), rB1
            #endif

         #endif
         #if KB > 159
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7536(pA), rA2
               prefetcht0 7536+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7536(pB), rB2
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
               movaps 7552(pB), rB0
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
               movaps 7568(pB), rB1
            #endif

         #endif
         #if KB > 160
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7584(pA), rA2
               prefetcht0 7584(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7584(pB), rB2
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
               movaps 7600(pB), rB0
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
               movaps 7616(pB), rB1
            #endif

         #endif
         #if KB > 161
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7632(pA), rA2
               prefetcht0 7632+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7632(pB), rB2
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
               movaps 7648(pB), rB0
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
               movaps 7664(pB), rB1
            #endif

         #endif
         #if KB > 162
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7680(pA), rA2
               prefetcht0 7680(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7680(pB), rB2
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
               movaps 7696(pB), rB0
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
               movaps 7712(pB), rB1
            #endif

         #endif
         #if KB > 163
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7728(pA), rA2
               prefetcht0 7728+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7728(pB), rB2
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
               movaps 7744(pB), rB0
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
               movaps 7760(pB), rB1
            #endif

         #endif
         #if KB > 164
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7776(pA), rA2
               prefetcht0 7776(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7776(pB), rB2
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
               movaps 7792(pB), rB0
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
               movaps 7808(pB), rB1
            #endif

         #endif
         #if KB > 165
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7824(pA), rA2
               prefetcht0 7824+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7824(pB), rB2
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
               movaps 7840(pB), rB0
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
               movaps 7856(pB), rB1
            #endif

         #endif
         #if KB > 166
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7872(pA), rA2
               prefetcht0 7872(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7872(pB), rB2
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
               movaps 7888(pB), rB0
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
               movaps 7904(pB), rB1
            #endif

         #endif
         #if KB > 167
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7920(pA), rA2
               prefetcht0 7920+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7920(pB), rB2
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
               movaps 7936(pB), rB0
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
               movaps 7952(pB), rB1
            #endif

         #endif
         #if KB > 168
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 7968(pA), rA2
               prefetcht0 7968(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 7968(pB), rB2
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
               movaps 7984(pB), rB0
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
               movaps 8000(pB), rB1
            #endif

         #endif
         #if KB > 169
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8016(pA), rA2
               prefetcht0 8016+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8016(pB), rB2
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
               movaps 8032(pB), rB0
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
               movaps 8048(pB), rB1
            #endif

         #endif
         #if KB > 170
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8064(pA), rA2
               prefetcht0 8064(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8064(pB), rB2
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
               movaps 8080(pB), rB0
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
               movaps 8096(pB), rB1
            #endif

         #endif
         #if KB > 171
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8112(pA), rA2
               prefetcht0 8112+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8112(pB), rB2
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
               movaps 8128(pB), rB0
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
               movaps 8144(pB), rB1
            #endif

         #endif
         #if KB > 172
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8160(pA), rA2
               prefetcht0 8160(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8160(pB), rB2
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
               movaps 8176(pB), rB0
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
               movaps 8192(pB), rB1
            #endif

         #endif
         #if KB > 173
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8208(pA), rA2
               prefetcht0 8208+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8208(pB), rB2
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
               movaps 8224(pB), rB0
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
               movaps 8240(pB), rB1
            #endif

         #endif
         #if KB > 174
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8256(pA), rA2
               prefetcht0 8256(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8256(pB), rB2
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
               movaps 8272(pB), rB0
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
               movaps 8288(pB), rB1
            #endif

         #endif
         #if KB > 175
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8304(pA), rA2
               prefetcht0 8304+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8304(pB), rB2
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
               movaps 8320(pB), rB0
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
               movaps 8336(pB), rB1
            #endif

         #endif
         #if KB > 176
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8352(pA), rA2
               prefetcht0 8352(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8352(pB), rB2
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
               movaps 8368(pB), rB0
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
               movaps 8384(pB), rB1
            #endif

         #endif
         #if KB > 177
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8400(pA), rA2
               prefetcht0 8400+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8400(pB), rB2
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
               movaps 8416(pB), rB0
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
               movaps 8432(pB), rB1
            #endif

         #endif
         #if KB > 178
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8448(pA), rA2
               prefetcht0 8448(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8448(pB), rB2
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
               movaps 8464(pB), rB0
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
               movaps 8480(pB), rB1
            #endif

         #endif
         #if KB > 179
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8496(pA), rA2
               prefetcht0 8496+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8496(pB), rB2
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
               movaps 8512(pB), rB0
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
               movaps 8528(pB), rB1
            #endif

         #endif
         #if KB > 180
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8544(pA), rA2
               prefetcht0 8544(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8544(pB), rB2
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
               movaps 8560(pB), rB0
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
               movaps 8576(pB), rB1
            #endif

         #endif
         #if KB > 181
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8592(pA), rA2
               prefetcht0 8592+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8592(pB), rB2
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
               movaps 8608(pB), rB0
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
               movaps 8624(pB), rB1
            #endif

         #endif
         #if KB > 182
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8640(pA), rA2
               prefetcht0 8640(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8640(pB), rB2
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
               movaps 8656(pB), rB0
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
               movaps 8672(pB), rB1
            #endif

         #endif
         #if KB > 183
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8688(pA), rA2
               prefetcht0 8688+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8688(pB), rB2
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
               movaps 8704(pB), rB0
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
               movaps 8720(pB), rB1
            #endif

         #endif
         #if KB > 184
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8736(pA), rA2
               prefetcht0 8736(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8736(pB), rB2
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
               movaps 8752(pB), rB0
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
               movaps 8768(pB), rB1
            #endif

         #endif
         #if KB > 185
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8784(pA), rA2
               prefetcht0 8784+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8784(pB), rB2
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
               movaps 8800(pB), rB0
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
               movaps 8816(pB), rB1
            #endif

         #endif
         #if KB > 186
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8832(pA), rA2
               prefetcht0 8832(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8832(pB), rB2
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
               movaps 8848(pB), rB0
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
               movaps 8864(pB), rB1
            #endif

         #endif
         #if KB > 187
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8880(pA), rA2
               prefetcht0 8880+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8880(pB), rB2
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
               movaps 8896(pB), rB0
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
               movaps 8912(pB), rB1
            #endif

         #endif
         #if KB > 188
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8928(pA), rA2
               prefetcht0 8928(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8928(pB), rB2
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
               movaps 8944(pB), rB0
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
               movaps 8960(pB), rB1
            #endif

         #endif
         #if KB > 189
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 8976(pA), rA2
               prefetcht0 8976+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 8976(pB), rB2
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
               movaps 8992(pB), rB0
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
               movaps 9008(pB), rB1
            #endif

         #endif
         #if KB > 190
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9024(pA), rA2
               prefetcht0 9024(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9024(pB), rB2
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
               movaps 9040(pB), rB0
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
               movaps 9056(pB), rB1
            #endif

         #endif
         #if KB > 191
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9072(pA), rA2
               prefetcht0 9072+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9072(pB), rB2
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
               movaps 9088(pB), rB0
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
               movaps 9104(pB), rB1
            #endif

         #endif
         #if KB > 192
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9120(pA), rA2
               prefetcht0 9120(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9120(pB), rB2
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
               movaps 9136(pB), rB0
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
               movaps 9152(pB), rB1
            #endif

         #endif
         #if KB > 193
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9168(pA), rA2
               prefetcht0 9168+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9168(pB), rB2
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
               movaps 9184(pB), rB0
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
               movaps 9200(pB), rB1
            #endif

         #endif
         #if KB > 194
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9216(pA), rA2
               prefetcht0 9216(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9216(pB), rB2
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
               movaps 9232(pB), rB0
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
               movaps 9248(pB), rB1
            #endif

         #endif
         #if KB > 195
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9264(pA), rA2
               prefetcht0 9264+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9264(pB), rB2
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
               movaps 9280(pB), rB0
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
               movaps 9296(pB), rB1
            #endif

         #endif
         #if KB > 196
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9312(pA), rA2
               prefetcht0 9312(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9312(pB), rB2
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
               movaps 9328(pB), rB0
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
               movaps 9344(pB), rB1
            #endif

         #endif
         #if KB > 197
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9360(pA), rA2
               prefetcht0 9360+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9360(pB), rB2
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
               movaps 9376(pB), rB0
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
               movaps 9392(pB), rB1
            #endif

         #endif
         #if KB > 198
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9408(pA), rA2
               prefetcht0 9408(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9408(pB), rB2
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
               movaps 9424(pB), rB0
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
               movaps 9440(pB), rB1
            #endif

         #endif
         #if KB > 199
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9456(pA), rA2
               prefetcht0 9456+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9456(pB), rB2
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
               movaps 9472(pB), rB0
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
               movaps 9488(pB), rB1
            #endif

         #endif
         #if KB > 200
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9504(pA), rA2
               prefetcht0 9504(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9504(pB), rB2
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
               movaps 9520(pB), rB0
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
               movaps 9536(pB), rB1
            #endif

         #endif
         #if KB > 201
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9552(pA), rA2
               prefetcht0 9552+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9552(pB), rB2
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
               movaps 9568(pB), rB0
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
               movaps 9584(pB), rB1
            #endif

         #endif
         #if KB > 202
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9600(pA), rA2
               prefetcht0 9600(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9600(pB), rB2
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
               movaps 9616(pB), rB0
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
               movaps 9632(pB), rB1
            #endif

         #endif
         #if KB > 203
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9648(pA), rA2
               prefetcht0 9648+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9648(pB), rB2
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
               movaps 9664(pB), rB0
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
               movaps 9680(pB), rB1
            #endif

         #endif
         #if KB > 204
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9696(pA), rA2
               prefetcht0 9696(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9696(pB), rB2
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
               movaps 9712(pB), rB0
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
               movaps 9728(pB), rB1
            #endif

         #endif
         #if KB > 205
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9744(pA), rA2
               prefetcht0 9744+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9744(pB), rB2
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
               movaps 9760(pB), rB0
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
               movaps 9776(pB), rB1
            #endif

         #endif
         #if KB > 206
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9792(pA), rA2
               prefetcht0 9792(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9792(pB), rB2
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
               movaps 9808(pB), rB0
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
               movaps 9824(pB), rB1
            #endif

         #endif
         #if KB > 207
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9840(pA), rA2
               prefetcht0 9840+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9840(pB), rB2
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
               movaps 9856(pB), rB0
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
               movaps 9872(pB), rB1
            #endif

         #endif
         #if KB > 208
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9888(pA), rA2
               prefetcht0 9888(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9888(pB), rB2
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
               movaps 9904(pB), rB0
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
               movaps 9920(pB), rB1
            #endif

         #endif
         #if KB > 209
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9936(pA), rA2
               prefetcht0 9936+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9936(pB), rB2
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
               movaps 9952(pB), rB0
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
               movaps 9968(pB), rB1
            #endif

         #endif
         #if KB > 210
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 9984(pA), rA2
               prefetcht0 9984(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 9984(pB), rB2
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
               movaps 10000(pB), rB0
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
               movaps 10016(pB), rB1
            #endif

         #endif
         #if KB > 211
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10032(pA), rA2
               prefetcht0 10032+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10032(pB), rB2
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
               movaps 10048(pB), rB0
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
               movaps 10064(pB), rB1
            #endif

         #endif
         #if KB > 212
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10080(pA), rA2
               prefetcht0 10080(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10080(pB), rB2
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
               movaps 10096(pB), rB0
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
               movaps 10112(pB), rB1
            #endif

         #endif
         #if KB > 213
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10128(pA), rA2
               prefetcht0 10128+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10128(pB), rB2
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
               movaps 10144(pB), rB0
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
               movaps 10160(pB), rB1
            #endif

         #endif
         #if KB > 214
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10176(pA), rA2
               prefetcht0 10176(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10176(pB), rB2
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
               movaps 10192(pB), rB0
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
               movaps 10208(pB), rB1
            #endif

         #endif
         #if KB > 215
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10224(pA), rA2
               prefetcht0 10224+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10224(pB), rB2
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
               movaps 10240(pB), rB0
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
               movaps 10256(pB), rB1
            #endif

         #endif
         #if KB > 216
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10272(pA), rA2
               prefetcht0 10272(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10272(pB), rB2
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
               movaps 10288(pB), rB0
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
               movaps 10304(pB), rB1
            #endif

         #endif
         #if KB > 217
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10320(pA), rA2
               prefetcht0 10320+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10320(pB), rB2
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
               movaps 10336(pB), rB0
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
               movaps 10352(pB), rB1
            #endif

         #endif
         #if KB > 218
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10368(pA), rA2
               prefetcht0 10368(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10368(pB), rB2
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
               movaps 10384(pB), rB0
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
               movaps 10400(pB), rB1
            #endif

         #endif
         #if KB > 219
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10416(pA), rA2
               prefetcht0 10416+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10416(pB), rB2
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
               movaps 10432(pB), rB0
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
               movaps 10448(pB), rB1
            #endif

         #endif
         #if KB > 220
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10464(pA), rA2
               prefetcht0 10464(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10464(pB), rB2
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
               movaps 10480(pB), rB0
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
               movaps 10496(pB), rB1
            #endif

         #endif
         #if KB > 221
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10512(pA), rA2
               prefetcht0 10512+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10512(pB), rB2
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
               movaps 10528(pB), rB0
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
               movaps 10544(pB), rB1
            #endif

         #endif
         #if KB > 222
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10560(pA), rA2
               prefetcht0 10560(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10560(pB), rB2
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
               movaps 10576(pB), rB0
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
               movaps 10592(pB), rB1
            #endif

         #endif
         #if KB > 223
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10608(pA), rA2
               prefetcht0 10608+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10608(pB), rB2
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
               movaps 10624(pB), rB0
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
               movaps 10640(pB), rB1
            #endif

         #endif
         #if KB > 224
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10656(pA), rA2
               prefetcht0 10656(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10656(pB), rB2
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
               movaps 10672(pB), rB0
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
               movaps 10688(pB), rB1
            #endif

         #endif
         #if KB > 225
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10704(pA), rA2
               prefetcht0 10704+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10704(pB), rB2
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
               movaps 10720(pB), rB0
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
               movaps 10736(pB), rB1
            #endif

         #endif
         #if KB > 226
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10752(pA), rA2
               prefetcht0 10752(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10752(pB), rB2
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
               movaps 10768(pB), rB0
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
               movaps 10784(pB), rB1
            #endif

         #endif
         #if KB > 227
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10800(pA), rA2
               prefetcht0 10800+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10800(pB), rB2
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
               movaps 10816(pB), rB0
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
               movaps 10832(pB), rB1
            #endif

         #endif
         #if KB > 228
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10848(pA), rA2
               prefetcht0 10848(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10848(pB), rB2
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
               movaps 10864(pB), rB0
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
               movaps 10880(pB), rB1
            #endif

         #endif
         #if KB > 229
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10896(pA), rA2
               prefetcht0 10896+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10896(pB), rB2
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
               movaps 10912(pB), rB0
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
               movaps 10928(pB), rB1
            #endif

         #endif
         #if KB > 230
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10944(pA), rA2
               prefetcht0 10944(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10944(pB), rB2
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
               movaps 10960(pB), rB0
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
               movaps 10976(pB), rB1
            #endif

         #endif
         #if KB > 231
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 10992(pA), rA2
               prefetcht0 10992+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 10992(pB), rB2
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
               movaps 11008(pB), rB0
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
               movaps 11024(pB), rB1
            #endif

         #endif
         #if KB > 232
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11040(pA), rA2
               prefetcht0 11040(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11040(pB), rB2
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
               movaps 11056(pB), rB0
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
               movaps 11072(pB), rB1
            #endif

         #endif
         #if KB > 233
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11088(pA), rA2
               prefetcht0 11088+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11088(pB), rB2
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
               movaps 11104(pB), rB0
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
               movaps 11120(pB), rB1
            #endif

         #endif
         #if KB > 234
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11136(pA), rA2
               prefetcht0 11136(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11136(pB), rB2
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
               movaps 11152(pB), rB0
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
               movaps 11168(pB), rB1
            #endif

         #endif
         #if KB > 235
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11184(pA), rA2
               prefetcht0 11184+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11184(pB), rB2
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
               movaps 11200(pB), rB0
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
               movaps 11216(pB), rB1
            #endif

         #endif
         #if KB > 236
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11232(pA), rA2
               prefetcht0 11232(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11232(pB), rB2
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
               movaps 11248(pB), rB0
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
               movaps 11264(pB), rB1
            #endif

         #endif
         #if KB > 237
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11280(pA), rA2
               prefetcht0 11280+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11280(pB), rB2
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
               movaps 11296(pB), rB0
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
               movaps 11312(pB), rB1
            #endif

         #endif
         #if KB > 238
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11328(pA), rA2
               prefetcht0 11328(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11328(pB), rB2
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
               movaps 11344(pB), rB0
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
               movaps 11360(pB), rB1
            #endif

         #endif
         #if KB > 239
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11376(pA), rA2
               prefetcht0 11376+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11376(pB), rB2
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
               movaps 11392(pB), rB0
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
               movaps 11408(pB), rB1
            #endif

         #endif
         #if KB > 240
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11424(pA), rA2
               prefetcht0 11424(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11424(pB), rB2
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
               movaps 11440(pB), rB0
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
               movaps 11456(pB), rB1
            #endif

         #endif
         #if KB > 241
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11472(pA), rA2
               prefetcht0 11472+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11472(pB), rB2
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
               movaps 11488(pB), rB0
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
               movaps 11504(pB), rB1
            #endif

         #endif
         #if KB > 242
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11520(pA), rA2
               prefetcht0 11520(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11520(pB), rB2
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
               movaps 11536(pB), rB0
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
               movaps 11552(pB), rB1
            #endif

         #endif
         #if KB > 243
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11568(pA), rA2
               prefetcht0 11568+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11568(pB), rB2
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
               movaps 11584(pB), rB0
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
               movaps 11600(pB), rB1
            #endif

         #endif
         #if KB > 244
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11616(pA), rA2
               prefetcht0 11616(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11616(pB), rB2
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
               movaps 11632(pB), rB0
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
               movaps 11648(pB), rB1
            #endif

         #endif
         #if KB > 245
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11664(pA), rA2
               prefetcht0 11664+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11664(pB), rB2
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
               movaps 11680(pB), rB0
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
               movaps 11696(pB), rB1
            #endif

         #endif
         #if KB > 246
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11712(pA), rA2
               prefetcht0 11712(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11712(pB), rB2
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
               movaps 11728(pB), rB0
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
               movaps 11744(pB), rB1
            #endif

         #endif
         #if KB > 247
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11760(pA), rA2
               prefetcht0 11760+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11760(pB), rB2
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
               movaps 11776(pB), rB0
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
               movaps 11792(pB), rB1
            #endif

         #endif
         #if KB > 248
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11808(pA), rA2
               prefetcht0 11808(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11808(pB), rB2
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
               movaps 11824(pB), rB0
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
               movaps 11840(pB), rB1
            #endif

         #endif
         #if KB > 249
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11856(pA), rA2
               prefetcht0 11856+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11856(pB), rB2
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
               movaps 11872(pB), rB0
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
               movaps 11888(pB), rB1
            #endif

         #endif
         #if KB > 250
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11904(pA), rA2
               prefetcht0 11904(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11904(pB), rB2
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
               movaps 11920(pB), rB0
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
               movaps 11936(pB), rB1
            #endif

         #endif
         #if KB > 251
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 11952(pA), rA2
               prefetcht0 11952+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 11952(pB), rB2
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
               movaps 11968(pB), rB0
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
               movaps 11984(pB), rB1
            #endif

         #endif
         #if KB > 252
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 12000(pA), rA2
               prefetcht0 12000(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 12000(pB), rB2
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
               movaps 12016(pB), rB0
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
               movaps 12032(pB), rB1
            #endif

         #endif
         #if KB > 253
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 12048(pA), rA2
               prefetcht0 12048+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 12048(pB), rB2
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
               movaps 12064(pB), rB0
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
               movaps 12080(pB), rB1
            #endif

         #endif
         #if KB > 254
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 12096(pA), rA2
               prefetcht0 12096(pA,incAm)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 12096(pB), rB2
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
               movaps 12112(pB), rB0
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
               movaps 12128(pB), rB1
            #endif

         #endif
         #if KB > 255
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            movaps 12144(pA), rA2
               prefetcht0 12144+PFBDIST(pB)

            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10

            movaps 12144(pB), rB2
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
               movaps 12160(pB), rB0
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
               movaps 12176(pB), rB1
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
