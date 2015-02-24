/*
 *             Automatically Tuned Linear Algebra Software v3.11.31
 * Copyright (C) 2012 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "atlas_asm.h"
#ifndef KB
   #define KB 0
#endif
/*
 * innermost (K-) loop items get priority on 1st 7 regs
 */
#define pA      %rcx
#define pB      %rdi
#define KK0     %rax
#define KK      %rdx
#define incA    %rsi
#define incB    %rbx
#define pA0     %r11
/*
 * Second (N-) loop items get next level of priority on good regs
 */
#define pC      %rbp
#define pfA     %r8
#define pfB     %r9
#define incPF   %r10
#define nnu     %r12
/*
 * Outer- (M-) loop variables assigned to any regs
 */
#define nmu     %r13
#define pB0     %r14
#define nnu0    %r15
/*
 * floating point registers
 */
#define m0   %xmm0
#define rA0  %xmm1
#define rA1  %xmm2
#define rA2  %xmm3
#define rB0  %xmm4
#define rB1  %xmm5
#define rB2  %xmm6
#define rC00 %xmm7
#define rC10 %xmm8
#define rC20 %xmm9
#define rC01 %xmm10
#define rC11 %xmm11
#define rC21 %xmm12
#define rC02 %xmm13
#define rC12 %xmm14
#define rC22 %xmm15
/*
                    rdi      rsi    rdx        rcx         r8        r9
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                  8(%rsp)    16(%rsp)     24(%rsp)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */

#define FSIZE 6*8
#ifndef prefA
   #define prefA prefetcht0
#endif
#ifndef prefB
   #define prefB prefetcht2
#endif
#ifndef prefC
   #ifdef ATL_3DNow
      #define prefC prefetchw
   #else
      #define prefC prefetcht0
   #endif
#endif
#ifdef BETAN1
   #define BETCOP subpd
#else
   #define BETCOP addpd
#endif
#define movapd movaps
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
/*
 * Save callee-saved iregs
 */
   sub $FSIZE, %rsp
   movq    %rbp, 0(%rsp)
   movq    %rbx, 8(%rsp)
   movq    %r12, 16(%rsp)
   movq    %r13, 24(%rsp)
   movq    %r14, 32(%rsp)
   movq    %r15, 40(%rsp)
/*
 * Load paramaters
 */
   mov %rdi, nmu
   mov %rsi, nnu
   mov %r8, pB
   mov %r9, pC
   mov nnu, nnu0
   movq FSIZE+8(%rsp), pfB      /* pfB = pAn */
   movq FSIZE+16(%rsp), pfA     /* pf = pBn */
   cmp pfA, pB
   CMOVE pfB, pfA
   CMOVEq FSIZE+24(%rsp), pfB
   mov $6*3*8, incPF           /* incPF = mu*nu*sizeof */
/*
 * Extend range of 1-byte offsets  by starting at -128
 */
   sub $-128, pA
   sub $-128, pB
   sub $-128, pC
   sub $-128, pfA
   sub $-128, pfB
   movq pB, pB0
   mov pA, pA0
   mov KK, KK0
   mov $192, incA
   mov $96, incB
            movddup -128(pB), rB0
            movddup -120(pB), rB1
   nop ; nop ; nop ; nop ; nop
   .local MNLOOP
   MNLOOP:
/*
      .local NLOOP
      NLOOP:
*/
/*
 *       Peel first iteration of K loop to initialize rCx
 */
         movapd -128(pA), rC02
         movapd rC02, rC00
         mulpd rB0, rC00
         movapd -112(pA), rC12
         movapd rC12, rC10
         mulpd rB0, rC10
         movapd -96(pA), rC22
         movapd rC22, rC20
         mulpd rB0, rC20

         movapd rC02, rC01
         mulpd rB1, rC01
            movddup -112(pB), rB2
         movapd rC12, rC11
         mulpd rB1, rC11
            prefC -128(pC)
         movapd rC22, rC21
         mulpd rB1, rC21

         mulpd rB2, rC02
         #if KB > 1
            movddup -104(pB), rB0
         #endif
         mulpd rB2, rC12
         #if KB > 1
            movddup -96(pB), rB1
         #endif
         mulpd rB2, rC22

         movapd -80(pA), rA0
         movapd rA0, m0
         mulpd rB0, m0
         addpd m0, rC00
         movapd -64(pA), rA1
         movapd rA1, m0
         mulpd rB0, m0
         addpd m0, rC10
         movapd -48(pA), rA2
         mulpd rA2, rB0
         addpd rB0, rC20

            prefC 64(pC)
         movapd rA0, m0
         mulpd rB1, m0
         addpd m0, rC01
            prefA -128(pfA)
         movapd rA1, m0
         mulpd rB1, m0
         addpd m0, rC11
         mulpd rA2, rB1
         addpd rB1, rC21

         movddup -88(pB), rB2
         mulpd rB2, rA0
         addpd rA0, rC02
            movddup -80(pB), rB0
         mulpd rB2, rA1
         addpd rA1, rC12
            movddup -72(pB), rB1
         mulpd rB2, rA2
         addpd rA2, rC22

         movapd -32(pA), rA0
         movapd rA0, m0
         mulpd rB0, m0
         addpd m0, rC00
         movapd -16(pA), rA1
         movapd rA1, m0
         mulpd rB0, m0
         addpd m0, rC10
         movapd (pA), rA2
         mulpd rA2, rB0
         addpd rB0, rC20

            prefA -64(pfA)
         movapd rA0, m0
         mulpd rB1, m0
         addpd m0, rC01
            prefA (pfA)
         movapd rA1, m0
         mulpd rB1, m0
         addpd m0, rC11
         mulpd rA2, rB1
         addpd rB1, rC21

         movddup -64(pB), rB2
         mulpd rB2, rA0
         addpd rA0, rC02
            movddup -56(pB), rB0
         mulpd rB2, rA1
         addpd rA1, rC12
            movddup -48(pB), rB1
         mulpd rB2, rA2
         addpd rA2, rC22

         movapd 16(pA), rA0
         movapd rA0, m0
         mulpd rB0, m0
         addpd m0, rC00
         movapd 32(pA), rA1
         movapd rA1, m0
         mulpd rB0, m0
         addpd m0, rC10
         movapd 48(pA), rA2
         mulpd rA2, rB0
         addpd rB0, rC20

            prefB -128(pfB)
         movapd rA0, m0
         mulpd rB1, m0
         addpd m0, rC01
            prefB -64(pfB)
         movapd rA1, m0
         mulpd rB1, m0
         addpd m0, rC11
            prefB (pfB)
         mulpd rA2, rB1
         addpd rB1, rC21

         movddup -40(pB), rB2
         add incB, pB
         mulpd rB2, rA0
         addpd rA0, rC02
            movddup -128(pB), rB0
         mulpd rB2, rA1
         addpd rA1, rC12
            movddup -120(pB), rB1
         mulpd rB2, rA2
            add incA, pA
         addpd rA2, rC22
         sub $4, KK
         jz KLOOPDONE

         KLOOP:
            movapd -128(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -112(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -96(pA), rA2
            mulpd rA2, rB0
            addpd rB0, rC20

            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup -112(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
               movddup -104(pB), rB0
            mulpd rB2, rA1
            addpd rA1, rC12
               movddup -96(pB), rB1
            mulpd rB2, rA2
            addpd rA2, rC22

            movapd -80(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -64(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -48(pA), rA2
            mulpd rA2, rB0
            addpd rB0, rC20

            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup -88(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
               movddup -80(pB), rB0
            mulpd rB2, rA1
            addpd rA1, rC12
               movddup -72(pB), rB1
            mulpd rB2, rA2
            addpd rA2, rC22

            movapd -32(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -16(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd (pA), rA2
            mulpd rA2, rB0
            addpd rB0, rC20

            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup -64(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
               movddup -56(pB), rB0
            mulpd rB2, rA1
            addpd rA1, rC12
               movddup -48(pB), rB1
            mulpd rB2, rA2
            addpd rA2, rC22

            movapd 16(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 32(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 48(pA), rA2
            mulpd rA2, rB0
            addpd rB0, rC20

            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
               add incA, pA
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup -40(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
               movddup -32(pB), rB0
            mulpd rB2, rA1
            addpd rA1, rC12
               movddup -24(pB), rB1
            mulpd rB2, rA2
            add incB, pB
            addpd rA2, rC22
         sub $4, KK
         jnz KLOOP
         KLOOPDONE:
         add incPF, pfA
         add incPF, pfB
/*
 *       Write answer back out to C
 */
         #ifdef BETA0
            movapd rC00, -128(pC)
            movapd rC10, -112(pC)
            movapd rC20, -96(pC)
            movapd rC01, -80(pC)
            movapd rC11, -64(pC)
            movapd rC21, -48(pC)
            movapd rC02, -32(pC)
            movapd rC12, -16(pC)
            movapd rC22, (pC)
/*
 *          Add running sum in rCx with original C, then store back out
 */
         #else
            BETCOP -128(pC), rC00
            movapd rC00, -128(pC)
            BETCOP -112(pC), rC10
            movapd rC10, -112(pC)
            BETCOP -96(pC), rC20
            movapd rC20, -96(pC)
            BETCOP -80(pC), rC01
            movapd rC01, -80(pC)
            BETCOP -64(pC), rC11
            movapd rC11, -64(pC)
            BETCOP -48(pC), rC21
            movapd rC21, -48(pC)
            BETCOP -32(pC), rC02
            movapd rC02, -32(pC)
            BETCOP -16(pC), rC12
            movapd rC12, -16(pC)
            BETCOP (pC), rC22
            movapd rC22, (pC)
         #endif
         add $6*3*8, pC        /* pC += MU*NU*sizeof */
         movddup -128(pB), rB0
      sub $1, nnu
         movddup -120(pB), rB1
         mov KK0, KK
         mov pA0, pA
      jnz MNLOOP

         movddup -128(pB0), rB0
      mov nnu0, nnu
         movddup -120(pB0), rB1
      mov pB0, pB
      lea (KK, KK, 4), pA      /* pA = 5*K */
      add KK, pA               /* pA = 6*K = MU*K */
      shl $3, pA               /* pA = MU*K*sizeof */
      add pA0, pA              /* pA = A + MU*K*sizeof */
      mov pA, pA0
   sub $1, nmu
   jnz MNLOOP
/* DONE: */
   movq    (%rsp), %rbp
   movq    8(%rsp), %rbx
   movq    16(%rsp), %r12
   movq    24(%rsp), %r13
   movq    32(%rsp), %r14
   movq    40(%rsp), %r15
   add $FSIZE, %rsp
   ret
