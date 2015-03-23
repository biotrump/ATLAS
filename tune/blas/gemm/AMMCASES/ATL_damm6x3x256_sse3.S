/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
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
#define i256    %rax
#define i768    %rdx   /* 3 * 256 */
#define i1280   %rsi   /* 5 * 256 */
#define i1792   %rbx   /* 7 * 256 */
#define i2304   %r11
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
   mov $256, i256
   lea (i256, i256,2), i768
   lea (i256, i256,4), i1280
   lea (i256, i768,2), i1792
   lea (i256, i256,8), i2304
            movddup -128(pB), rB0
            movddup -120(pB), rB1
   ALIGN8
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

/*
 *       Fully unrolled K-loop
 */
         #if KB > 1
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

               prefC (pC)
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
            #if KB > 2
            movddup -80(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 2
            movddup -72(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 2
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

               prefA -128(pfA)
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
            #if KB > 3
            movddup -56(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 3
            movddup -48(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 3
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

               prefA -64(pfA)
            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup -40(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 4
            movddup -32(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 4
            movddup -24(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 4
            movapd 64(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 80(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 96(pA), rA2
            mulpd rA2, rB0
            addpd rB0, rC20

               prefA (pfA)
            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup -16(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 5
            movddup -8(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 5
            movddup (pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 5
            movapd 112(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -128(pA,i256), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -112(pA,i256), rA2
            mulpd rA2, rB0
            addpd rB0, rC20

               prefB -128(pfB)
            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup 8(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 6
            movddup 16(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 6
            movddup 24(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 6
            movapd -96(pA,i256), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -80(pA,i256), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -64(pA,i256), rA2
            mulpd rA2, rB0
            addpd rB0, rC20

               prefB -64(pfB)
            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup 32(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 7
            movddup 40(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 7
            movddup 48(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 7
            movapd -48(pA,i256), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -32(pA,i256), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -16(pA,i256), rA2
            mulpd rA2, rB0
            addpd rB0, rC20

               prefB (pfB)
            movapd rA0, m0
            mulpd rB1, m0
            addpd m0, rC01
            movapd rA1, m0
            mulpd rB1, m0
            addpd m0, rC11
            mulpd rA2, rB1
            addpd rB1, rC21

            movddup 56(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 8
            movddup 64(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 8
            movddup 72(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 8
            movapd (pA,i256), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 16(pA,i256), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 32(pA,i256), rA2
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

            movddup 80(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 9
            movddup 88(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 9
            movddup 96(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 9
            movapd 48(pA,i256), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 64(pA,i256), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 80(pA,i256), rA2
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

            movddup 104(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 10
            movddup 112(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 10
            movddup 120(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 10
            movapd 96(pA,i256), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 112(pA,i256), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -128(pA,i256,2), rA2
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

            movddup -128(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 11
            movddup -120(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 11
            movddup -112(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 11
            movapd -112(pA,i256,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -96(pA,i256,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -80(pA,i256,2), rA2
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

            movddup -104(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 12
            movddup -96(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 12
            movddup -88(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 12
            movapd -64(pA,i256,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -48(pA,i256,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -32(pA,i256,2), rA2
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

            movddup -80(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 13
            movddup -72(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 13
            movddup -64(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 13
            movapd -16(pA,i256,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd (pA,i256,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 16(pA,i256,2), rA2
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

            movddup -56(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 14
            movddup -48(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 14
            movddup -40(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 14
            movapd 32(pA,i256,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 48(pA,i256,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 64(pA,i256,2), rA2
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

            movddup -32(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 15
            movddup -24(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 15
            movddup -16(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 15
            movapd 80(pA,i256,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 96(pA,i256,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 112(pA,i256,2), rA2
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

            movddup -8(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 16
            movddup (pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 16
            movddup 8(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 16
            movapd -128(pA,i768), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -112(pA,i768), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -96(pA,i768), rA2
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

            movddup 16(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 17
            movddup 24(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 17
            movddup 32(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 17
            movapd -80(pA,i768), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -64(pA,i768), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -48(pA,i768), rA2
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

            movddup 40(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 18
            movddup 48(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 18
            movddup 56(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 18
            movapd -32(pA,i768), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -16(pA,i768), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd (pA,i768), rA2
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

            movddup 64(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 19
            movddup 72(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 19
            movddup 80(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 19
            movapd 16(pA,i768), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 32(pA,i768), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 48(pA,i768), rA2
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

            movddup 88(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 20
            movddup 96(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 20
            movddup 104(pB,i256), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 20
            movapd 64(pA,i768), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 80(pA,i768), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 96(pA,i768), rA2
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

            movddup 112(pB,i256), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 21
            movddup 120(pB,i256), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 21
            movddup -128(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 21
            movapd 112(pA,i768), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -128(pA,i256,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -112(pA,i256,4), rA2
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

            movddup -120(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 22
            movddup -112(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 22
            movddup -104(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 22
            movapd -96(pA,i256,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -80(pA,i256,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -64(pA,i256,4), rA2
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

            movddup -96(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 23
            movddup -88(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 23
            movddup -80(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 23
            movapd -48(pA,i256,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -32(pA,i256,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -16(pA,i256,4), rA2
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

            movddup -72(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 24
            movddup -64(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 24
            movddup -56(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 24
            movapd (pA,i256,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 16(pA,i256,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 32(pA,i256,4), rA2
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

            movddup -48(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 25
            movddup -40(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 25
            movddup -32(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 25
            movapd 48(pA,i256,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 64(pA,i256,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 80(pA,i256,4), rA2
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

            movddup -24(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 26
            movddup -16(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 26
            movddup -8(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 26
            movapd 96(pA,i256,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 112(pA,i256,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -128(pA,i1280), rA2
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

            movddup (pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 27
            movddup 8(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 27
            movddup 16(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 27
            movapd -112(pA,i1280), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -96(pA,i1280), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -80(pA,i1280), rA2
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

            movddup 24(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 28
            movddup 32(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 28
            movddup 40(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 28
            movapd -64(pA,i1280), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -48(pA,i1280), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -32(pA,i1280), rA2
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

            movddup 48(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 29
            movddup 56(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 29
            movddup 64(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 29
            movapd -16(pA,i1280), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd (pA,i1280), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 16(pA,i1280), rA2
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

            movddup 72(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 30
            movddup 80(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 30
            movddup 88(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 30
            movapd 32(pA,i1280), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 48(pA,i1280), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 64(pA,i1280), rA2
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

            movddup 96(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 31
            movddup 104(pB,i256,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 31
            movddup 112(pB,i256,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 31
            movapd 80(pA,i1280), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 96(pA,i1280), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 112(pA,i1280), rA2
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

            movddup 120(pB,i256,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 32
            movddup -128(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 32
            movddup -120(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 32
            movapd -128(pA,i768,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -112(pA,i768,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -96(pA,i768,2), rA2
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

            movddup -112(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 33
            movddup -104(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 33
            movddup -96(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 33
            movapd -80(pA,i768,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -64(pA,i768,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -48(pA,i768,2), rA2
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

            movddup -88(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 34
            movddup -80(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 34
            movddup -72(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 34
            movapd -32(pA,i768,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -16(pA,i768,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd (pA,i768,2), rA2
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

            movddup -64(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 35
            movddup -56(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 35
            movddup -48(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 35
            movapd 16(pA,i768,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 32(pA,i768,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 48(pA,i768,2), rA2
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

            movddup -40(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 36
            movddup -32(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 36
            movddup -24(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 36
            movapd 64(pA,i768,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 80(pA,i768,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 96(pA,i768,2), rA2
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

            movddup -16(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 37
            movddup -8(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 37
            movddup (pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 37
            movapd 112(pA,i768,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -128(pA,i1792), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -112(pA,i1792), rA2
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

            movddup 8(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 38
            movddup 16(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 38
            movddup 24(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 38
            movapd -96(pA,i1792), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -80(pA,i1792), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -64(pA,i1792), rA2
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

            movddup 32(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 39
            movddup 40(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 39
            movddup 48(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 39
            movapd -48(pA,i1792), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -32(pA,i1792), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -16(pA,i1792), rA2
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

            movddup 56(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 40
            movddup 64(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 40
            movddup 72(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 40
            movapd (pA,i1792), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 16(pA,i1792), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 32(pA,i1792), rA2
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

            movddup 80(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 41
            movddup 88(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 41
            movddup 96(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 41
            movapd 48(pA,i1792), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 64(pA,i1792), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 80(pA,i1792), rA2
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

            movddup 104(pB,i768), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 42
            movddup 112(pB,i768), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 42
            movddup 120(pB,i768), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 42
            movapd 96(pA,i1792), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 112(pA,i1792), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -128(pA,i256,8), rA2
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

            movddup -128(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 43
            movddup -120(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 43
            movddup -112(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 43
            movapd -112(pA,i256,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -96(pA,i256,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -80(pA,i256,8), rA2
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

            movddup -104(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 44
            movddup -96(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 44
            movddup -88(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 44
            movapd -64(pA,i256,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -48(pA,i256,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -32(pA,i256,8), rA2
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

            movddup -80(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 45
            movddup -72(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 45
            movddup -64(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 45
            movapd -16(pA,i256,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd (pA,i256,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 16(pA,i256,8), rA2
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

            movddup -56(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 46
            movddup -48(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 46
            movddup -40(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 46
            movapd 32(pA,i256,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 48(pA,i256,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 64(pA,i256,8), rA2
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

            movddup -32(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 47
            movddup -24(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 47
            movddup -16(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 47
            movapd 80(pA,i256,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 96(pA,i256,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 112(pA,i256,8), rA2
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

            movddup -8(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 48
            movddup (pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 48
            movddup 8(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 48
            movapd -128(pA,i2304), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -112(pA,i2304), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -96(pA,i2304), rA2
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

            movddup 16(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 49
            movddup 24(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 49
            movddup 32(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 49
            movapd -80(pA,i2304), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -64(pA,i2304), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -48(pA,i2304), rA2
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

            movddup 40(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 50
            movddup 48(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 50
            movddup 56(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 50
            movapd -32(pA,i2304), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -16(pA,i2304), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd (pA,i2304), rA2
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

            movddup 64(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 51
            movddup 72(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 51
            movddup 80(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 51
            movapd 16(pA,i2304), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 32(pA,i2304), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 48(pA,i2304), rA2
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

            movddup 88(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 52
            movddup 96(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 52
            movddup 104(pB,i256,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 52
            movapd 64(pA,i2304), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 80(pA,i2304), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 96(pA,i2304), rA2
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

            movddup 112(pB,i256,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 53
            movddup 120(pB,i256,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 53
            movddup -128(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 53
            movapd 112(pA,i2304), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -128(pA,i1280,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -112(pA,i1280,2), rA2
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

            movddup -120(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 54
            movddup -112(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 54
            movddup -104(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 54
            movapd -96(pA,i1280,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -80(pA,i1280,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -64(pA,i1280,2), rA2
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

            movddup -96(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 55
            movddup -88(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 55
            movddup -80(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 55
            movapd -48(pA,i1280,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -32(pA,i1280,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -16(pA,i1280,2), rA2
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

            movddup -72(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 56
            movddup -64(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 56
            movddup -56(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 56
            movapd (pA,i1280,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 16(pA,i1280,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 32(pA,i1280,2), rA2
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

            movddup -48(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 57
            movddup -40(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 57
            movddup -32(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 57
            movapd 48(pA,i1280,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 64(pA,i1280,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 80(pA,i1280,2), rA2
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

            movddup -24(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 58
            movddup -16(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 58
            movddup -8(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 58
            movapd 96(pA,i1280,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 112(pA,i1280,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 2688(pA), rA2
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

            movddup (pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 59
            movddup 8(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 59
            movddup 16(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 59
            movapd 2704(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 2720(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 2736(pA), rA2
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

            movddup 24(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 60
            movddup 32(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 60
            movddup 40(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 60
            movapd 2752(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 2768(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 2784(pA), rA2
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

            movddup 48(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 61
            movddup 56(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 61
            movddup 64(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 61
            movapd 2800(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 2816(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 2832(pA), rA2
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

            movddup 72(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 62
            movddup 80(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 62
            movddup 88(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 62
            movapd 2848(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 2864(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 2880(pA), rA2
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

            movddup 96(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 63
            movddup 104(pB,i1280), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 63
            movddup 112(pB,i1280), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 63
            movapd 2896(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 2912(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 2928(pA), rA2
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

            movddup 120(pB,i1280), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 64
            movddup -128(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 64
            movddup -120(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 64
            movapd -128(pA,i768,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -112(pA,i768,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -96(pA,i768,4), rA2
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

            movddup -112(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 65
            movddup -104(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 65
            movddup -96(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 65
            movapd -80(pA,i768,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -64(pA,i768,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -48(pA,i768,4), rA2
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

            movddup -88(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 66
            movddup -80(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 66
            movddup -72(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 66
            movapd -32(pA,i768,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -16(pA,i768,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd (pA,i768,4), rA2
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

            movddup -64(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 67
            movddup -56(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 67
            movddup -48(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 67
            movapd 16(pA,i768,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 32(pA,i768,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 48(pA,i768,4), rA2
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

            movddup -40(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 68
            movddup -32(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 68
            movddup -24(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 68
            movapd 64(pA,i768,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 80(pA,i768,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 96(pA,i768,4), rA2
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

            movddup -16(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 69
            movddup -8(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 69
            movddup (pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 69
            movapd 112(pA,i768,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3200(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3216(pA), rA2
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

            movddup 8(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 70
            movddup 16(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 70
            movddup 24(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 70
            movapd 3232(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3248(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3264(pA), rA2
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

            movddup 32(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 71
            movddup 40(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 71
            movddup 48(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 71
            movapd 3280(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3296(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3312(pA), rA2
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

            movddup 56(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 72
            movddup 64(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 72
            movddup 72(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 72
            movapd 3328(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3344(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3360(pA), rA2
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

            movddup 80(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 73
            movddup 88(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 73
            movddup 96(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 73
            movapd 3376(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3392(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3408(pA), rA2
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

            movddup 104(pB,i768,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 74
            movddup 112(pB,i768,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 74
            movddup 120(pB,i768,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 74
            movapd 3424(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3440(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -128(pA,i1792,2), rA2
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

            movddup -128(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 75
            movddup -120(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 75
            movddup -112(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 75
            movapd -112(pA,i1792,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -96(pA,i1792,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -80(pA,i1792,2), rA2
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

            movddup -104(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 76
            movddup -96(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 76
            movddup -88(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 76
            movapd -64(pA,i1792,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -48(pA,i1792,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -32(pA,i1792,2), rA2
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

            movddup -80(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 77
            movddup -72(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 77
            movddup -64(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 77
            movapd -16(pA,i1792,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd (pA,i1792,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 16(pA,i1792,2), rA2
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

            movddup -56(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 78
            movddup -48(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 78
            movddup -40(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 78
            movapd 32(pA,i1792,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 48(pA,i1792,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 64(pA,i1792,2), rA2
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

            movddup -32(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 79
            movddup -24(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 79
            movddup -16(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 79
            movapd 80(pA,i1792,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 96(pA,i1792,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 112(pA,i1792,2), rA2
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

            movddup -8(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 80
            movddup (pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 80
            movddup 8(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 80
            movapd 3712(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3728(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3744(pA), rA2
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

            movddup 16(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 81
            movddup 24(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 81
            movddup 32(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 81
            movapd 3760(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3776(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3792(pA), rA2
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

            movddup 40(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 82
            movddup 48(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 82
            movddup 56(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 82
            movapd 3808(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3824(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3840(pA), rA2
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

            movddup 64(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 83
            movddup 72(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 83
            movddup 80(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 83
            movapd 3856(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3872(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3888(pA), rA2
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

            movddup 88(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 84
            movddup 96(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 84
            movddup 104(pB,i1792), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 84
            movapd 3904(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3920(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3936(pA), rA2
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

            movddup 112(pB,i1792), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 85
            movddup 120(pB,i1792), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 85
            movddup -128(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 85
            movapd 3952(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 3968(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 3984(pA), rA2
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

            movddup -120(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 86
            movddup -112(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 86
            movddup -104(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 86
            movapd 4000(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4016(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4032(pA), rA2
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

            movddup -96(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 87
            movddup -88(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 87
            movddup -80(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 87
            movapd 4048(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4064(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4080(pA), rA2
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

            movddup -72(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 88
            movddup -64(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 88
            movddup -56(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 88
            movapd 4096(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4112(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4128(pA), rA2
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

            movddup -48(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 89
            movddup -40(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 89
            movddup -32(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 89
            movapd 4144(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4160(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4176(pA), rA2
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

            movddup -24(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 90
            movddup -16(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 90
            movddup -8(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 90
            movapd 4192(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4208(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4224(pA), rA2
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

            movddup (pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 91
            movddup 8(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 91
            movddup 16(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 91
            movapd 4240(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4256(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4272(pA), rA2
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

            movddup 24(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 92
            movddup 32(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 92
            movddup 40(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 92
            movapd 4288(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4304(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4320(pA), rA2
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

            movddup 48(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 93
            movddup 56(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 93
            movddup 64(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 93
            movapd 4336(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4352(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4368(pA), rA2
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

            movddup 72(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 94
            movddup 80(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 94
            movddup 88(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 94
            movapd 4384(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4400(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4416(pA), rA2
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

            movddup 96(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 95
            movddup 104(pB,i256,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 95
            movddup 112(pB,i256,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 95
            movapd 4432(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4448(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4464(pA), rA2
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

            movddup 120(pB,i256,8), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 96
            movddup -128(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 96
            movddup -120(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 96
            movapd -128(pA,i2304,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -112(pA,i2304,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -96(pA,i2304,2), rA2
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

            movddup -112(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 97
            movddup -104(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 97
            movddup -96(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 97
            movapd -80(pA,i2304,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -64(pA,i2304,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -48(pA,i2304,2), rA2
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

            movddup -88(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 98
            movddup -80(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 98
            movddup -72(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 98
            movapd -32(pA,i2304,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -16(pA,i2304,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd (pA,i2304,2), rA2
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

            movddup -64(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 99
            movddup -56(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 99
            movddup -48(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 99
            movapd 16(pA,i2304,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 32(pA,i2304,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 48(pA,i2304,2), rA2
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

            movddup -40(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 100
            movddup -32(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 100
            movddup -24(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 100
            movapd 64(pA,i2304,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 80(pA,i2304,2), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 96(pA,i2304,2), rA2
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

            movddup -16(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 101
            movddup -8(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 101
            movddup (pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 101
            movapd 112(pA,i2304,2), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4736(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4752(pA), rA2
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

            movddup 8(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 102
            movddup 16(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 102
            movddup 24(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 102
            movapd 4768(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4784(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4800(pA), rA2
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

            movddup 32(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 103
            movddup 40(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 103
            movddup 48(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 103
            movapd 4816(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4832(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4848(pA), rA2
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

            movddup 56(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 104
            movddup 64(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 104
            movddup 72(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 104
            movapd 4864(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4880(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4896(pA), rA2
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

            movddup 80(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 105
            movddup 88(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 105
            movddup 96(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 105
            movapd 4912(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4928(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 4944(pA), rA2
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

            movddup 104(pB,i2304), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 106
            movddup 112(pB,i2304), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 106
            movddup 120(pB,i2304), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 106
            movapd 4960(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 4976(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -128(pA,i1280,4), rA2
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

            movddup -128(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 107
            movddup -120(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 107
            movddup -112(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 107
            movapd -112(pA,i1280,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -96(pA,i1280,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -80(pA,i1280,4), rA2
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

            movddup -104(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 108
            movddup -96(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 108
            movddup -88(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 108
            movapd -64(pA,i1280,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -48(pA,i1280,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -32(pA,i1280,4), rA2
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

            movddup -80(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 109
            movddup -72(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 109
            movddup -64(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 109
            movapd -16(pA,i1280,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd (pA,i1280,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 16(pA,i1280,4), rA2
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

            movddup -56(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 110
            movddup -48(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 110
            movddup -40(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 110
            movapd 32(pA,i1280,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 48(pA,i1280,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 64(pA,i1280,4), rA2
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

            movddup -32(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 111
            movddup -24(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 111
            movddup -16(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 111
            movapd 80(pA,i1280,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 96(pA,i1280,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 112(pA,i1280,4), rA2
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

            movddup -8(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 112
            movddup (pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 112
            movddup 8(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 112
            movapd 5248(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5264(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5280(pA), rA2
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

            movddup 16(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 113
            movddup 24(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 113
            movddup 32(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 113
            movapd 5296(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5312(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5328(pA), rA2
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

            movddup 40(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 114
            movddup 48(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 114
            movddup 56(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 114
            movapd 5344(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5360(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5376(pA), rA2
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

            movddup 64(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 115
            movddup 72(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 115
            movddup 80(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 115
            movapd 5392(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5408(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5424(pA), rA2
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

            movddup 88(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 116
            movddup 96(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 116
            movddup 104(pB,i1280,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 116
            movapd 5440(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5456(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5472(pA), rA2
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

            movddup 112(pB,i1280,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 117
            movddup 120(pB,i1280,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 117
            movddup 2688(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 117
            movapd 5488(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5504(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5520(pA), rA2
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

            movddup 2696(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 118
            movddup 2704(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 118
            movddup 2712(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 118
            movapd 5536(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5552(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5568(pA), rA2
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

            movddup 2720(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 119
            movddup 2728(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 119
            movddup 2736(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 119
            movapd 5584(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5600(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5616(pA), rA2
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

            movddup 2744(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 120
            movddup 2752(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 120
            movddup 2760(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 120
            movapd 5632(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5648(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5664(pA), rA2
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

            movddup 2768(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 121
            movddup 2776(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 121
            movddup 2784(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 121
            movapd 5680(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5696(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5712(pA), rA2
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

            movddup 2792(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 122
            movddup 2800(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 122
            movddup 2808(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 122
            movapd 5728(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5744(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5760(pA), rA2
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

            movddup 2816(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 123
            movddup 2824(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 123
            movddup 2832(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 123
            movapd 5776(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5792(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5808(pA), rA2
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

            movddup 2840(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 124
            movddup 2848(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 124
            movddup 2856(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 124
            movapd 5824(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5840(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5856(pA), rA2
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

            movddup 2864(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 125
            movddup 2872(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 125
            movddup 2880(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 125
            movapd 5872(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5888(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5904(pA), rA2
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

            movddup 2888(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 126
            movddup 2896(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 126
            movddup 2904(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 126
            movapd 5920(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5936(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 5952(pA), rA2
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

            movddup 2912(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 127
            movddup 2920(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 127
            movddup 2928(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 127
            movapd 5968(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 5984(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6000(pA), rA2
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

            movddup 2936(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 128
            movddup -128(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 128
            movddup -120(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 128
            movapd -128(pA,i768,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -112(pA,i768,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -96(pA,i768,8), rA2
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

            movddup -112(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 129
            movddup -104(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 129
            movddup -96(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 129
            movapd -80(pA,i768,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -64(pA,i768,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -48(pA,i768,8), rA2
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

            movddup -88(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 130
            movddup -80(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 130
            movddup -72(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 130
            movapd -32(pA,i768,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -16(pA,i768,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd (pA,i768,8), rA2
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

            movddup -64(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 131
            movddup -56(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 131
            movddup -48(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 131
            movapd 16(pA,i768,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 32(pA,i768,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 48(pA,i768,8), rA2
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

            movddup -40(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 132
            movddup -32(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 132
            movddup -24(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 132
            movapd 64(pA,i768,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 80(pA,i768,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 96(pA,i768,8), rA2
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

            movddup -16(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 133
            movddup -8(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 133
            movddup (pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 133
            movapd 112(pA,i768,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6272(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6288(pA), rA2
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

            movddup 8(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 134
            movddup 16(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 134
            movddup 24(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 134
            movapd 6304(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6320(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6336(pA), rA2
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

            movddup 32(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 135
            movddup 40(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 135
            movddup 48(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 135
            movapd 6352(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6368(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6384(pA), rA2
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

            movddup 56(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 136
            movddup 64(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 136
            movddup 72(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 136
            movapd 6400(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6416(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6432(pA), rA2
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

            movddup 80(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 137
            movddup 88(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 137
            movddup 96(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 137
            movapd 6448(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6464(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6480(pA), rA2
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

            movddup 104(pB,i768,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 138
            movddup 112(pB,i768,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 138
            movddup 120(pB,i768,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 138
            movapd 6496(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6512(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6528(pA), rA2
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

            movddup 3200(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 139
            movddup 3208(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 139
            movddup 3216(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 139
            movapd 6544(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6560(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6576(pA), rA2
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

            movddup 3224(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 140
            movddup 3232(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 140
            movddup 3240(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 140
            movapd 6592(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6608(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6624(pA), rA2
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

            movddup 3248(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 141
            movddup 3256(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 141
            movddup 3264(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 141
            movapd 6640(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6656(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6672(pA), rA2
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

            movddup 3272(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 142
            movddup 3280(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 142
            movddup 3288(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 142
            movapd 6688(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6704(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6720(pA), rA2
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

            movddup 3296(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 143
            movddup 3304(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 143
            movddup 3312(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 143
            movapd 6736(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6752(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6768(pA), rA2
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

            movddup 3320(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 144
            movddup 3328(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 144
            movddup 3336(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 144
            movapd 6784(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6800(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6816(pA), rA2
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

            movddup 3344(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 145
            movddup 3352(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 145
            movddup 3360(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 145
            movapd 6832(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6848(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6864(pA), rA2
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

            movddup 3368(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 146
            movddup 3376(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 146
            movddup 3384(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 146
            movapd 6880(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6896(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6912(pA), rA2
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

            movddup 3392(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 147
            movddup 3400(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 147
            movddup 3408(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 147
            movapd 6928(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6944(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 6960(pA), rA2
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

            movddup 3416(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 148
            movddup 3424(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 148
            movddup 3432(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 148
            movapd 6976(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 6992(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7008(pA), rA2
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

            movddup 3440(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 149
            movddup 3448(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 149
            movddup -128(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 149
            movapd 7024(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -128(pA,i1792,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -112(pA,i1792,4), rA2
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

            movddup -120(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 150
            movddup -112(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 150
            movddup -104(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 150
            movapd -96(pA,i1792,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -80(pA,i1792,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -64(pA,i1792,4), rA2
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

            movddup -96(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 151
            movddup -88(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 151
            movddup -80(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 151
            movapd -48(pA,i1792,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -32(pA,i1792,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -16(pA,i1792,4), rA2
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

            movddup -72(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 152
            movddup -64(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 152
            movddup -56(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 152
            movapd (pA,i1792,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 16(pA,i1792,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 32(pA,i1792,4), rA2
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

            movddup -48(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 153
            movddup -40(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 153
            movddup -32(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 153
            movapd 48(pA,i1792,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 64(pA,i1792,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 80(pA,i1792,4), rA2
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

            movddup -24(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 154
            movddup -16(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 154
            movddup -8(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 154
            movapd 96(pA,i1792,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 112(pA,i1792,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7296(pA), rA2
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

            movddup (pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 155
            movddup 8(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 155
            movddup 16(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 155
            movapd 7312(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7328(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7344(pA), rA2
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

            movddup 24(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 156
            movddup 32(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 156
            movddup 40(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 156
            movapd 7360(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7376(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7392(pA), rA2
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

            movddup 48(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 157
            movddup 56(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 157
            movddup 64(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 157
            movapd 7408(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7424(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7440(pA), rA2
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

            movddup 72(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 158
            movddup 80(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 158
            movddup 88(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 158
            movapd 7456(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7472(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7488(pA), rA2
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

            movddup 96(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 159
            movddup 104(pB,i1792,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 159
            movddup 112(pB,i1792,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 159
            movapd 7504(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7520(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7536(pA), rA2
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

            movddup 120(pB,i1792,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 160
            movddup 3712(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 160
            movddup 3720(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 160
            movapd 7552(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7568(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7584(pA), rA2
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

            movddup 3728(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 161
            movddup 3736(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 161
            movddup 3744(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 161
            movapd 7600(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7616(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7632(pA), rA2
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

            movddup 3752(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 162
            movddup 3760(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 162
            movddup 3768(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 162
            movapd 7648(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7664(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7680(pA), rA2
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

            movddup 3776(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 163
            movddup 3784(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 163
            movddup 3792(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 163
            movapd 7696(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7712(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7728(pA), rA2
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

            movddup 3800(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 164
            movddup 3808(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 164
            movddup 3816(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 164
            movapd 7744(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7760(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7776(pA), rA2
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

            movddup 3824(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 165
            movddup 3832(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 165
            movddup 3840(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 165
            movapd 7792(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7808(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7824(pA), rA2
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

            movddup 3848(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 166
            movddup 3856(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 166
            movddup 3864(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 166
            movapd 7840(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7856(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7872(pA), rA2
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

            movddup 3872(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 167
            movddup 3880(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 167
            movddup 3888(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 167
            movapd 7888(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7904(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7920(pA), rA2
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

            movddup 3896(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 168
            movddup 3904(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 168
            movddup 3912(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 168
            movapd 7936(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 7952(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 7968(pA), rA2
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

            movddup 3920(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 169
            movddup 3928(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 169
            movddup 3936(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 169
            movapd 7984(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8000(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8016(pA), rA2
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

            movddup 3944(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 170
            movddup 3952(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 170
            movddup 3960(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 170
            movapd 8032(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8048(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8064(pA), rA2
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

            movddup 3968(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 171
            movddup 3976(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 171
            movddup 3984(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 171
            movapd 8080(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8096(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8112(pA), rA2
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

            movddup 3992(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 172
            movddup 4000(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 172
            movddup 4008(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 172
            movapd 8128(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8144(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8160(pA), rA2
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

            movddup 4016(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 173
            movddup 4024(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 173
            movddup 4032(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 173
            movapd 8176(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8192(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8208(pA), rA2
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

            movddup 4040(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 174
            movddup 4048(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 174
            movddup 4056(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 174
            movapd 8224(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8240(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8256(pA), rA2
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

            movddup 4064(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 175
            movddup 4072(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 175
            movddup 4080(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 175
            movapd 8272(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8288(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8304(pA), rA2
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

            movddup 4088(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 176
            movddup 4096(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 176
            movddup 4104(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 176
            movapd 8320(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8336(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8352(pA), rA2
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

            movddup 4112(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 177
            movddup 4120(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 177
            movddup 4128(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 177
            movapd 8368(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8384(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8400(pA), rA2
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

            movddup 4136(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 178
            movddup 4144(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 178
            movddup 4152(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 178
            movapd 8416(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8432(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8448(pA), rA2
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

            movddup 4160(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 179
            movddup 4168(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 179
            movddup 4176(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 179
            movapd 8464(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8480(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8496(pA), rA2
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

            movddup 4184(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 180
            movddup 4192(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 180
            movddup 4200(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 180
            movapd 8512(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8528(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8544(pA), rA2
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

            movddup 4208(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 181
            movddup 4216(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 181
            movddup 4224(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 181
            movapd 8560(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8576(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8592(pA), rA2
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

            movddup 4232(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 182
            movddup 4240(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 182
            movddup 4248(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 182
            movapd 8608(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8624(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8640(pA), rA2
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

            movddup 4256(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 183
            movddup 4264(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 183
            movddup 4272(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 183
            movapd 8656(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8672(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8688(pA), rA2
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

            movddup 4280(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 184
            movddup 4288(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 184
            movddup 4296(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 184
            movapd 8704(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8720(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8736(pA), rA2
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

            movddup 4304(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 185
            movddup 4312(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 185
            movddup 4320(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 185
            movapd 8752(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8768(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8784(pA), rA2
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

            movddup 4328(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 186
            movddup 4336(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 186
            movddup 4344(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 186
            movapd 8800(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8816(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8832(pA), rA2
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

            movddup 4352(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 187
            movddup 4360(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 187
            movddup 4368(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 187
            movapd 8848(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8864(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8880(pA), rA2
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

            movddup 4376(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 188
            movddup 4384(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 188
            movddup 4392(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 188
            movapd 8896(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8912(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8928(pA), rA2
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

            movddup 4400(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 189
            movddup 4408(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 189
            movddup 4416(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 189
            movapd 8944(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 8960(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 8976(pA), rA2
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

            movddup 4424(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 190
            movddup 4432(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 190
            movddup 4440(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 190
            movapd 8992(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9008(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9024(pA), rA2
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

            movddup 4448(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 191
            movddup 4456(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 191
            movddup 4464(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 191
            movapd 9040(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9056(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9072(pA), rA2
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

            movddup 4472(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 192
            movddup -128(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 192
            movddup -120(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 192
            movapd -128(pA,i2304,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -112(pA,i2304,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -96(pA,i2304,4), rA2
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

            movddup -112(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 193
            movddup -104(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 193
            movddup -96(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 193
            movapd -80(pA,i2304,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -64(pA,i2304,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -48(pA,i2304,4), rA2
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

            movddup -88(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 194
            movddup -80(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 194
            movddup -72(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 194
            movapd -32(pA,i2304,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -16(pA,i2304,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd (pA,i2304,4), rA2
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

            movddup -64(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 195
            movddup -56(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 195
            movddup -48(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 195
            movapd 16(pA,i2304,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 32(pA,i2304,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 48(pA,i2304,4), rA2
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

            movddup -40(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 196
            movddup -32(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 196
            movddup -24(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 196
            movapd 64(pA,i2304,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 80(pA,i2304,4), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 96(pA,i2304,4), rA2
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

            movddup -16(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 197
            movddup -8(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 197
            movddup (pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 197
            movapd 112(pA,i2304,4), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9344(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9360(pA), rA2
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

            movddup 8(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 198
            movddup 16(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 198
            movddup 24(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 198
            movapd 9376(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9392(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9408(pA), rA2
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

            movddup 32(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 199
            movddup 40(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 199
            movddup 48(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 199
            movapd 9424(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9440(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9456(pA), rA2
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

            movddup 56(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 200
            movddup 64(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 200
            movddup 72(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 200
            movapd 9472(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9488(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9504(pA), rA2
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

            movddup 80(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 201
            movddup 88(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 201
            movddup 96(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 201
            movapd 9520(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9536(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9552(pA), rA2
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

            movddup 104(pB,i2304,2), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 202
            movddup 112(pB,i2304,2), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 202
            movddup 120(pB,i2304,2), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 202
            movapd 9568(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9584(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9600(pA), rA2
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

            movddup 4736(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 203
            movddup 4744(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 203
            movddup 4752(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 203
            movapd 9616(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9632(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9648(pA), rA2
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

            movddup 4760(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 204
            movddup 4768(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 204
            movddup 4776(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 204
            movapd 9664(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9680(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9696(pA), rA2
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

            movddup 4784(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 205
            movddup 4792(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 205
            movddup 4800(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 205
            movapd 9712(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9728(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9744(pA), rA2
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

            movddup 4808(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 206
            movddup 4816(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 206
            movddup 4824(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 206
            movapd 9760(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9776(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9792(pA), rA2
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

            movddup 4832(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 207
            movddup 4840(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 207
            movddup 4848(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 207
            movapd 9808(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9824(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9840(pA), rA2
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

            movddup 4856(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 208
            movddup 4864(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 208
            movddup 4872(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 208
            movapd 9856(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9872(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9888(pA), rA2
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

            movddup 4880(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 209
            movddup 4888(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 209
            movddup 4896(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 209
            movapd 9904(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9920(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9936(pA), rA2
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

            movddup 4904(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 210
            movddup 4912(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 210
            movddup 4920(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 210
            movapd 9952(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 9968(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 9984(pA), rA2
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

            movddup 4928(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 211
            movddup 4936(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 211
            movddup 4944(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 211
            movapd 10000(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10016(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10032(pA), rA2
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

            movddup 4952(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 212
            movddup 4960(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 212
            movddup 4968(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 212
            movapd 10048(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10064(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10080(pA), rA2
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

            movddup 4976(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 213
            movddup 4984(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 213
            movddup -128(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 213
            movapd 10096(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -128(pA,i1280,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -112(pA,i1280,8), rA2
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

            movddup -120(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 214
            movddup -112(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 214
            movddup -104(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 214
            movapd -96(pA,i1280,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -80(pA,i1280,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -64(pA,i1280,8), rA2
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

            movddup -96(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 215
            movddup -88(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 215
            movddup -80(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 215
            movapd -48(pA,i1280,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd -32(pA,i1280,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd -16(pA,i1280,8), rA2
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

            movddup -72(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 216
            movddup -64(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 216
            movddup -56(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 216
            movapd (pA,i1280,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 16(pA,i1280,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 32(pA,i1280,8), rA2
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

            movddup -48(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 217
            movddup -40(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 217
            movddup -32(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 217
            movapd 48(pA,i1280,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 64(pA,i1280,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 80(pA,i1280,8), rA2
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

            movddup -24(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 218
            movddup -16(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 218
            movddup -8(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 218
            movapd 96(pA,i1280,8), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 112(pA,i1280,8), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10368(pA), rA2
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

            movddup (pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 219
            movddup 8(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 219
            movddup 16(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 219
            movapd 10384(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10400(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10416(pA), rA2
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

            movddup 24(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 220
            movddup 32(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 220
            movddup 40(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 220
            movapd 10432(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10448(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10464(pA), rA2
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

            movddup 48(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 221
            movddup 56(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 221
            movddup 64(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 221
            movapd 10480(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10496(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10512(pA), rA2
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

            movddup 72(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 222
            movddup 80(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 222
            movddup 88(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 222
            movapd 10528(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10544(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10560(pA), rA2
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

            movddup 96(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 223
            movddup 104(pB,i1280,4), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 223
            movddup 112(pB,i1280,4), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 223
            movapd 10576(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10592(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10608(pA), rA2
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

            movddup 120(pB,i1280,4), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 224
            movddup 5248(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 224
            movddup 5256(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 224
            movapd 10624(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10640(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10656(pA), rA2
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

            movddup 5264(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 225
            movddup 5272(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 225
            movddup 5280(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 225
            movapd 10672(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10688(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10704(pA), rA2
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

            movddup 5288(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 226
            movddup 5296(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 226
            movddup 5304(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 226
            movapd 10720(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10736(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10752(pA), rA2
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

            movddup 5312(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 227
            movddup 5320(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 227
            movddup 5328(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 227
            movapd 10768(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10784(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10800(pA), rA2
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

            movddup 5336(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 228
            movddup 5344(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 228
            movddup 5352(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 228
            movapd 10816(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10832(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10848(pA), rA2
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

            movddup 5360(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 229
            movddup 5368(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 229
            movddup 5376(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 229
            movapd 10864(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10880(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10896(pA), rA2
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

            movddup 5384(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 230
            movddup 5392(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 230
            movddup 5400(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 230
            movapd 10912(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10928(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10944(pA), rA2
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

            movddup 5408(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 231
            movddup 5416(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 231
            movddup 5424(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 231
            movapd 10960(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 10976(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 10992(pA), rA2
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

            movddup 5432(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 232
            movddup 5440(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 232
            movddup 5448(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 232
            movapd 11008(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11024(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11040(pA), rA2
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

            movddup 5456(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 233
            movddup 5464(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 233
            movddup 5472(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 233
            movapd 11056(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11072(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11088(pA), rA2
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

            movddup 5480(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 234
            movddup 5488(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 234
            movddup 5496(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 234
            movapd 11104(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11120(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11136(pA), rA2
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

            movddup 5504(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 235
            movddup 5512(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 235
            movddup 5520(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 235
            movapd 11152(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11168(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11184(pA), rA2
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

            movddup 5528(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 236
            movddup 5536(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 236
            movddup 5544(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 236
            movapd 11200(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11216(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11232(pA), rA2
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

            movddup 5552(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 237
            movddup 5560(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 237
            movddup 5568(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 237
            movapd 11248(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11264(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11280(pA), rA2
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

            movddup 5576(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 238
            movddup 5584(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 238
            movddup 5592(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 238
            movapd 11296(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11312(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11328(pA), rA2
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

            movddup 5600(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 239
            movddup 5608(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 239
            movddup 5616(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 239
            movapd 11344(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11360(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11376(pA), rA2
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

            movddup 5624(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 240
            movddup 5632(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 240
            movddup 5640(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 240
            movapd 11392(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11408(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11424(pA), rA2
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

            movddup 5648(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 241
            movddup 5656(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 241
            movddup 5664(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 241
            movapd 11440(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11456(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11472(pA), rA2
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

            movddup 5672(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 242
            movddup 5680(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 242
            movddup 5688(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 242
            movapd 11488(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11504(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11520(pA), rA2
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

            movddup 5696(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 243
            movddup 5704(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 243
            movddup 5712(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 243
            movapd 11536(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11552(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11568(pA), rA2
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

            movddup 5720(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 244
            movddup 5728(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 244
            movddup 5736(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 244
            movapd 11584(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11600(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11616(pA), rA2
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

            movddup 5744(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 245
            movddup 5752(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 245
            movddup 5760(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 245
            movapd 11632(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11648(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11664(pA), rA2
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

            movddup 5768(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 246
            movddup 5776(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 246
            movddup 5784(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 246
            movapd 11680(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11696(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11712(pA), rA2
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

            movddup 5792(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 247
            movddup 5800(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 247
            movddup 5808(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 247
            movapd 11728(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11744(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11760(pA), rA2
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

            movddup 5816(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 248
            movddup 5824(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 248
            movddup 5832(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 248
            movapd 11776(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11792(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11808(pA), rA2
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

            movddup 5840(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 249
            movddup 5848(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 249
            movddup 5856(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 249
            movapd 11824(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11840(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11856(pA), rA2
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

            movddup 5864(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 250
            movddup 5872(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 250
            movddup 5880(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 250
            movapd 11872(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11888(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11904(pA), rA2
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

            movddup 5888(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 251
            movddup 5896(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 251
            movddup 5904(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 251
            movapd 11920(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11936(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 11952(pA), rA2
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

            movddup 5912(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 252
            movddup 5920(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 252
            movddup 5928(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 252
            movapd 11968(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 11984(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 12000(pA), rA2
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

            movddup 5936(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 253
            movddup 5944(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 253
            movddup 5952(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 253
            movapd 12016(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 12032(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 12048(pA), rA2
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

            movddup 5960(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 254
            movddup 5968(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 254
            movddup 5976(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 254
            movapd 12064(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 12080(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 12096(pA), rA2
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

            movddup 5984(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 255
            movddup 5992(pB), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 255
            movddup 6000(pB), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
         #if KB > 255
            movapd 12112(pA), rA0
            movapd rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            movapd 12128(pA), rA1
            movapd rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            movapd 12144(pA), rA2
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

            movddup 6008(pB), rB2
            mulpd rB2, rA0
            addpd rA0, rC02
            #if KB > 256
            movddup -128(pB,i768,8), rB0
            #endif
            mulpd rB2, rA1
            addpd rA1, rC12
            #if KB > 256
            movddup -120(pB,i768,8), rB1
            #endif
            mulpd rB2, rA2
            addpd rA2, rC22
         #endif
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
         add $KB*3*8, pB        /* pB += K*NU*sizeof */
         add $6*3*8, pC        /* pC += MU*NU*sizeof */
         movddup -128(pB), rB0
      sub $1, nnu
         movddup -120(pB), rB1
      jnz MNLOOP

         movddup -128(pB0), rB0
      mov nnu0, nnu
         movddup -120(pB0), rB1
      mov pB0, pB
      add $KB*6*8, pA          /* pA += KB*MU*size */
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
