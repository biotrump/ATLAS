/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2013 R. Clint Whaley
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
#define movaps movaps
#define nmu     %rdi
#define nnu     %rsi
#define nnu0    %r10
#define pA      %rcx
#define pB      %rax
#define pC      %r9
#define pf      %rbp
#define pB0     %r12
#define incPF   %rbx
#define pfB     %rdx
#define incAm   %r11
#define pfC     %r13
#define r256    %r14

#define rm0     %xmm0
#define rA0     %xmm1
#define rA1     %xmm2
#define rB0     %xmm3
#define rB1     %xmm4
#define rB2     %xmm5
#define rB3     %xmm6
#define rb3     %xmm7
#define rC00    %xmm8
#define rC10    %xmm9
#define rC01    %xmm10
#define rC11    %xmm11
#define rC02    %xmm12
#define rC12    %xmm13
#define rC03    %xmm14
#define rC13    %xmm15
#ifndef pref
   #define pref prefetcht1
#endif
#ifndef prefB
   #define prefB prefetcht1
#endif
#ifndef prefC
   #ifdef ATL_3DNow
      #define prefC prefetchw
   #else
      #define prefC prefetcht0
   #endif
#endif
#ifdef BETAN1
   #define BETCOP subps
#else
   #define BETCOP addps
#endif
#define FSIZE 6*8
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
/*
 * Load paramaters
 */
   movq %r8, pB
   mov nnu, nnu0
   movq FSIZE+16(%rsp), pf      /* pf = pBn */
   movq FSIZE+8(%rsp), pfB      /* pfB = pAn */
   movq FSIZE+24(%rsp), pfC    /* pfC = pCn */
   mov $4*4*4, incPF
/*
 * Extend range of small operands by starting at -128
 */
   sub $-128, pA
   sub $-128, pB
   mov $KB*8*4, incAm           /* incAm = KB*MU*size */
   movq pB, pB0

   ALIGN8
   .local MNLOOP
   MNLOOP:
/*
 *       Peel first iteration of K-loop to handle init of C to 0
 */
         movaps -128(pB), rB3
         pshufd $0x00, rB3, rB0
         movaps -128(pA), rC00

         movaps rC00, rC01
         mulps rB0, rC00
         pshufd $0x55, rB3, rB1
         movaps -112(pA), rC10
         movaps rC10, rC11
         mulps rB0, rC10

         #if KB > 1
            movaps -112(pB), rb3
         #else
            pref (pf)
         #endif
         movaps rC01, rC02
         mulps rB1, rC01
         #if KB > 1
            movaps -96(pA), rA0
         #else
            pref 64(pf)
         #endif
         movaps rC11, rC12
         mulps rB1, rC11

         #if KB > 1
            movaps -80(pA), rA1
         #endif
         pshufd $0xAA, rB3, rB2
         movaps rC02, rC03
         mulps rB2, rC02
         shufps $0xFF, rB3, rB3
         movaps rC12, rC13
         mulps rB2, rC12
         #if KB > 1
            pshufd $0x00, rb3, rB0
         #else
            add incPF, pf
         #endif

         mulps rB3, rC03
         prefC (pC)
         mulps rB3, rC13
         prefC 64(pC)
/*
 *       2nd peeled K iteration
 */
         #if KB > 1
            #if KB > 2
               movaps -96(pB), rB3
            #endif
            pshufd $0x55, rb3, rB1
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            pshufd $0xAA, rb3, rB2
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            shufps $0xFF, rb3, rb3
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB > 2
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 2
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 2
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 2
               movaps -64(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 2
               movaps -48(pA), rA1
            #endif
         #endif
/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 2
            #if KB > 3
               movaps -80(pB), rb3
            #elif KB == 3
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 3
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 3
               pshufd $0x00, rb3, rB0
            #elif KB == 3
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 3
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 3
               pshufd $0x55, rb3, rB1
            #elif KB == 3
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 3
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 3
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 3
               movaps -32(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 3
               movaps -16(pA), rA1
            #endif
         #endif
         #if KB > 3
            #if KB > 4
               movaps -64(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 4
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 4
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 4
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 4
               movaps 0(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 4
               movaps 16(pA), rA1
            #endif
         #endif
         #if KB > 4
            #if KB > 5
               movaps -48(pB), rb3
            #elif KB == 5
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 5
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 5
               pshufd $0x00, rb3, rB0
            #elif KB == 5
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 5
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 5
               pshufd $0x55, rb3, rB1
            #elif KB == 5
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 5
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 5
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 5
               movaps 32(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 5
               movaps 48(pA), rA1
            #endif
         #endif
         #if KB > 5
            #if KB > 6
               movaps -32(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 6
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 6
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 6
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 6
               movaps 64(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 6
               movaps 80(pA), rA1
            #endif
         #endif
         #if KB > 6
            #if KB > 7
               movaps -16(pB), rb3
            #elif KB == 7
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 7
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 7
               pshufd $0x00, rb3, rB0
            #elif KB == 7
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 7
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 7
               pshufd $0x55, rb3, rB1
            #elif KB == 7
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 7
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 7
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 7
               movaps 96(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 7
               movaps 112(pA), rA1
            #endif
         #endif
         #if KB > 7
            #if KB > 8
               movaps 0(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 8
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 8
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 8
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 8
               movaps 128(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 8
               movaps 144(pA), rA1
            #endif
         #endif
         #if KB > 8
            #if KB > 9
               movaps 16(pB), rb3
            #elif KB == 9
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 9
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 9
               pshufd $0x00, rb3, rB0
            #elif KB == 9
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 9
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 9
               pshufd $0x55, rb3, rB1
            #elif KB == 9
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 9
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 9
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 9
               movaps 160(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 9
               movaps 176(pA), rA1
            #endif
         #endif
         #if KB > 9
            #if KB > 10
               movaps 32(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 10
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 10
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 10
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 10
               movaps 192(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 10
               movaps 208(pA), rA1
            #endif
         #endif
         #if KB > 10
            #if KB > 11
               movaps 48(pB), rb3
            #elif KB == 11
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 11
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 11
               pshufd $0x00, rb3, rB0
            #elif KB == 11
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 11
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 11
               pshufd $0x55, rb3, rB1
            #elif KB == 11
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 11
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 11
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 11
               movaps 224(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 11
               movaps 240(pA), rA1
            #endif
         #endif
         #if KB > 11
            #if KB > 12
               movaps 64(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 12
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 12
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 12
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 12
               movaps 256(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 12
               movaps 272(pA), rA1
            #endif
         #endif
         #if KB > 12
            #if KB > 13
               movaps 80(pB), rb3
            #elif KB == 13
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 13
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 13
               pshufd $0x00, rb3, rB0
            #elif KB == 13
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 13
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 13
               pshufd $0x55, rb3, rB1
            #elif KB == 13
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 13
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 13
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 13
               movaps 288(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 13
               movaps 304(pA), rA1
            #endif
         #endif
         #if KB > 13
            #if KB > 14
               movaps 96(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 14
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 14
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 14
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 14
               movaps 320(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 14
               movaps 336(pA), rA1
            #endif
         #endif
         #if KB > 14
            #if KB > 15
               movaps 112(pB), rb3
            #elif KB == 15
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 15
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 15
               pshufd $0x00, rb3, rB0
            #elif KB == 15
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 15
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 15
               pshufd $0x55, rb3, rB1
            #elif KB == 15
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 15
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 15
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 15
               movaps 352(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 15
               movaps 368(pA), rA1
            #endif
         #endif
         #if KB > 15
            #if KB > 16
               movaps 128(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 16
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 16
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 16
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 16
               movaps 384(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 16
               movaps 400(pA), rA1
            #endif
         #endif
         #if KB > 16
            #if KB > 17
               movaps 144(pB), rb3
            #elif KB == 17
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 17
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 17
               pshufd $0x00, rb3, rB0
            #elif KB == 17
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 17
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 17
               pshufd $0x55, rb3, rB1
            #elif KB == 17
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 17
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 17
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 17
               movaps 416(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 17
               movaps 432(pA), rA1
            #endif
         #endif
         #if KB > 17
            #if KB > 18
               movaps 160(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 18
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 18
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 18
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 18
               movaps 448(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 18
               movaps 464(pA), rA1
            #endif
         #endif
         #if KB > 18
            #if KB > 19
               movaps 176(pB), rb3
            #elif KB == 19
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 19
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 19
               pshufd $0x00, rb3, rB0
            #elif KB == 19
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 19
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 19
               pshufd $0x55, rb3, rB1
            #elif KB == 19
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 19
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 19
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 19
               movaps 480(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 19
               movaps 496(pA), rA1
            #endif
         #endif
         #if KB > 19
            #if KB > 20
               movaps 192(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 20
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 20
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 20
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 20
               movaps 512(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 20
               movaps 528(pA), rA1
            #endif
         #endif
         #if KB > 20
            #if KB > 21
               movaps 208(pB), rb3
            #elif KB == 21
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 21
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 21
               pshufd $0x00, rb3, rB0
            #elif KB == 21
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 21
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 21
               pshufd $0x55, rb3, rB1
            #elif KB == 21
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 21
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 21
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 21
               movaps 544(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 21
               movaps 560(pA), rA1
            #endif
         #endif
         #if KB > 21
            #if KB > 22
               movaps 224(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 22
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 22
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 22
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 22
               movaps 576(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 22
               movaps 592(pA), rA1
            #endif
         #endif
         #if KB > 22
            #if KB > 23
               movaps 240(pB), rb3
            #elif KB == 23
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 23
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 23
               pshufd $0x00, rb3, rB0
            #elif KB == 23
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 23
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 23
               pshufd $0x55, rb3, rB1
            #elif KB == 23
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 23
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 23
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 23
               movaps 608(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 23
               movaps 624(pA), rA1
            #endif
         #endif
         #if KB > 23
            #if KB > 24
               movaps 256(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 24
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 24
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 24
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 24
               movaps 640(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 24
               movaps 656(pA), rA1
            #endif
         #endif
         #if KB > 24
            #if KB > 25
               movaps 272(pB), rb3
            #elif KB == 25
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 25
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 25
               pshufd $0x00, rb3, rB0
            #elif KB == 25
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 25
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 25
               pshufd $0x55, rb3, rB1
            #elif KB == 25
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 25
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 25
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 25
               movaps 672(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 25
               movaps 688(pA), rA1
            #endif
         #endif
         #if KB > 25
            #if KB > 26
               movaps 288(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 26
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 26
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 26
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 26
               movaps 704(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 26
               movaps 720(pA), rA1
            #endif
         #endif
         #if KB > 26
            #if KB > 27
               movaps 304(pB), rb3
            #elif KB == 27
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 27
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 27
               pshufd $0x00, rb3, rB0
            #elif KB == 27
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 27
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 27
               pshufd $0x55, rb3, rB1
            #elif KB == 27
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 27
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 27
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 27
               movaps 736(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 27
               movaps 752(pA), rA1
            #endif
         #endif
         #if KB > 27
            #if KB > 28
               movaps 320(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 28
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 28
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 28
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 28
               movaps 768(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 28
               movaps 784(pA), rA1
            #endif
         #endif
         #if KB > 28
            #if KB > 29
               movaps 336(pB), rb3
            #elif KB == 29
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 29
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 29
               pshufd $0x00, rb3, rB0
            #elif KB == 29
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 29
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 29
               pshufd $0x55, rb3, rB1
            #elif KB == 29
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 29
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 29
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 29
               movaps 800(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 29
               movaps 816(pA), rA1
            #endif
         #endif
         #if KB > 29
            #if KB > 30
               movaps 352(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 30
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 30
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 30
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 30
               movaps 832(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 30
               movaps 848(pA), rA1
            #endif
         #endif
         #if KB > 30
            #if KB > 31
               movaps 368(pB), rb3
            #elif KB == 31
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 31
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 31
               pshufd $0x00, rb3, rB0
            #elif KB == 31
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 31
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 31
               pshufd $0x55, rb3, rB1
            #elif KB == 31
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 31
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 31
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 31
               movaps 864(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 31
               movaps 880(pA), rA1
            #endif
         #endif
         #if KB > 31
            #if KB > 32
               movaps 384(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 32
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 32
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 32
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 32
               movaps 896(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 32
               movaps 912(pA), rA1
            #endif
         #endif
         #if KB > 32
            #if KB > 33
               movaps 400(pB), rb3
            #elif KB == 33
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 33
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 33
               pshufd $0x00, rb3, rB0
            #elif KB == 33
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 33
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 33
               pshufd $0x55, rb3, rB1
            #elif KB == 33
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 33
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 33
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 33
               movaps 928(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 33
               movaps 944(pA), rA1
            #endif
         #endif
         #if KB > 33
            #if KB > 34
               movaps 416(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 34
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 34
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 34
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 34
               movaps 960(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 34
               movaps 976(pA), rA1
            #endif
         #endif
         #if KB > 34
            #if KB > 35
               movaps 432(pB), rb3
            #elif KB == 35
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 35
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 35
               pshufd $0x00, rb3, rB0
            #elif KB == 35
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 35
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 35
               pshufd $0x55, rb3, rB1
            #elif KB == 35
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 35
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 35
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 35
               movaps 992(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 35
               movaps 1008(pA), rA1
            #endif
         #endif
         #if KB > 35
            #if KB > 36
               movaps 448(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 36
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 36
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 36
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 36
               movaps 1024(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 36
               movaps 1040(pA), rA1
            #endif
         #endif
         #if KB > 36
            #if KB > 37
               movaps 464(pB), rb3
            #elif KB == 37
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 37
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 37
               pshufd $0x00, rb3, rB0
            #elif KB == 37
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 37
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 37
               pshufd $0x55, rb3, rB1
            #elif KB == 37
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 37
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 37
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 37
               movaps 1056(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 37
               movaps 1072(pA), rA1
            #endif
         #endif
         #if KB > 37
            #if KB > 38
               movaps 480(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 38
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 38
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 38
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 38
               movaps 1088(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 38
               movaps 1104(pA), rA1
            #endif
         #endif
         #if KB > 38
            #if KB > 39
               movaps 496(pB), rb3
            #elif KB == 39
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 39
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 39
               pshufd $0x00, rb3, rB0
            #elif KB == 39
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 39
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 39
               pshufd $0x55, rb3, rB1
            #elif KB == 39
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 39
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 39
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 39
               movaps 1120(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 39
               movaps 1136(pA), rA1
            #endif
         #endif
         #if KB > 39
            #if KB > 40
               movaps 512(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 40
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 40
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 40
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 40
               movaps 1152(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 40
               movaps 1168(pA), rA1
            #endif
         #endif
         #if KB > 40
            #if KB > 41
               movaps 528(pB), rb3
            #elif KB == 41
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 41
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 41
               pshufd $0x00, rb3, rB0
            #elif KB == 41
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 41
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 41
               pshufd $0x55, rb3, rB1
            #elif KB == 41
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 41
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 41
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 41
               movaps 1184(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 41
               movaps 1200(pA), rA1
            #endif
         #endif
         #if KB > 41
            #if KB > 42
               movaps 544(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 42
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 42
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 42
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 42
               movaps 1216(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 42
               movaps 1232(pA), rA1
            #endif
         #endif
         #if KB > 42
            #if KB > 43
               movaps 560(pB), rb3
            #elif KB == 43
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 43
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 43
               pshufd $0x00, rb3, rB0
            #elif KB == 43
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 43
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 43
               pshufd $0x55, rb3, rB1
            #elif KB == 43
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 43
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 43
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 43
               movaps 1248(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 43
               movaps 1264(pA), rA1
            #endif
         #endif
         #if KB > 43
            #if KB > 44
               movaps 576(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 44
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 44
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 44
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 44
               movaps 1280(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 44
               movaps 1296(pA), rA1
            #endif
         #endif
         #if KB > 44
            #if KB > 45
               movaps 592(pB), rb3
            #elif KB == 45
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 45
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 45
               pshufd $0x00, rb3, rB0
            #elif KB == 45
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 45
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 45
               pshufd $0x55, rb3, rB1
            #elif KB == 45
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 45
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 45
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 45
               movaps 1312(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 45
               movaps 1328(pA), rA1
            #endif
         #endif
         #if KB > 45
            #if KB > 46
               movaps 608(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 46
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 46
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 46
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 46
               movaps 1344(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 46
               movaps 1360(pA), rA1
            #endif
         #endif
         #if KB > 46
            #if KB > 47
               movaps 624(pB), rb3
            #elif KB == 47
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 47
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 47
               pshufd $0x00, rb3, rB0
            #elif KB == 47
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 47
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 47
               pshufd $0x55, rb3, rB1
            #elif KB == 47
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 47
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 47
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 47
               movaps 1376(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 47
               movaps 1392(pA), rA1
            #endif
         #endif
         #if KB > 47
            #if KB > 48
               movaps 640(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 48
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 48
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 48
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 48
               movaps 1408(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 48
               movaps 1424(pA), rA1
            #endif
         #endif
         #if KB > 48
            #if KB > 49
               movaps 656(pB), rb3
            #elif KB == 49
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 49
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 49
               pshufd $0x00, rb3, rB0
            #elif KB == 49
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 49
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 49
               pshufd $0x55, rb3, rB1
            #elif KB == 49
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 49
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 49
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 49
               movaps 1440(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 49
               movaps 1456(pA), rA1
            #endif
         #endif
         #if KB > 49
            #if KB > 50
               movaps 672(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 50
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 50
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 50
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 50
               movaps 1472(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 50
               movaps 1488(pA), rA1
            #endif
         #endif
         #if KB > 50
            #if KB > 51
               movaps 688(pB), rb3
            #elif KB == 51
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 51
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 51
               pshufd $0x00, rb3, rB0
            #elif KB == 51
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 51
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 51
               pshufd $0x55, rb3, rB1
            #elif KB == 51
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 51
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 51
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 51
               movaps 1504(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 51
               movaps 1520(pA), rA1
            #endif
         #endif
         #if KB > 51
            #if KB > 52
               movaps 704(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 52
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 52
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 52
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 52
               movaps 1536(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 52
               movaps 1552(pA), rA1
            #endif
         #endif
         #if KB > 52
            #if KB > 53
               movaps 720(pB), rb3
            #elif KB == 53
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 53
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 53
               pshufd $0x00, rb3, rB0
            #elif KB == 53
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 53
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 53
               pshufd $0x55, rb3, rB1
            #elif KB == 53
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 53
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 53
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 53
               movaps 1568(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 53
               movaps 1584(pA), rA1
            #endif
         #endif
         #if KB > 53
            #if KB > 54
               movaps 736(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 54
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 54
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 54
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 54
               movaps 1600(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 54
               movaps 1616(pA), rA1
            #endif
         #endif
         #if KB > 54
            #if KB > 55
               movaps 752(pB), rb3
            #elif KB == 55
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 55
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 55
               pshufd $0x00, rb3, rB0
            #elif KB == 55
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 55
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 55
               pshufd $0x55, rb3, rB1
            #elif KB == 55
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 55
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 55
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 55
               movaps 1632(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 55
               movaps 1648(pA), rA1
            #endif
         #endif
         #if KB > 55
            #if KB > 56
               movaps 768(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 56
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 56
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 56
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 56
               movaps 1664(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 56
               movaps 1680(pA), rA1
            #endif
         #endif
         #if KB > 56
            #if KB > 57
               movaps 784(pB), rb3
            #elif KB == 57
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 57
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 57
               pshufd $0x00, rb3, rB0
            #elif KB == 57
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 57
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 57
               pshufd $0x55, rb3, rB1
            #elif KB == 57
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 57
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 57
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 57
               movaps 1696(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 57
               movaps 1712(pA), rA1
            #endif
         #endif
         #if KB > 57
            #if KB > 58
               movaps 800(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 58
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 58
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 58
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 58
               movaps 1728(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 58
               movaps 1744(pA), rA1
            #endif
         #endif
         #if KB > 58
            #if KB > 59
               movaps 816(pB), rb3
            #elif KB == 59
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 59
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 59
               pshufd $0x00, rb3, rB0
            #elif KB == 59
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 59
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 59
               pshufd $0x55, rb3, rB1
            #elif KB == 59
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 59
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 59
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 59
               movaps 1760(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 59
               movaps 1776(pA), rA1
            #endif
         #endif
         #if KB > 59
            #if KB > 60
               movaps 832(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 60
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 60
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 60
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 60
               movaps 1792(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 60
               movaps 1808(pA), rA1
            #endif
         #endif
         #if KB > 60
            #if KB > 61
               movaps 848(pB), rb3
            #elif KB == 61
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 61
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 61
               pshufd $0x00, rb3, rB0
            #elif KB == 61
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 61
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 61
               pshufd $0x55, rb3, rB1
            #elif KB == 61
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 61
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 61
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 61
               movaps 1824(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 61
               movaps 1840(pA), rA1
            #endif
         #endif
         #if KB > 61
            #if KB > 62
               movaps 864(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 62
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 62
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 62
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 62
               movaps 1856(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 62
               movaps 1872(pA), rA1
            #endif
         #endif
         #if KB > 62
            #if KB > 63
               movaps 880(pB), rb3
            #elif KB == 63
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 63
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 63
               pshufd $0x00, rb3, rB0
            #elif KB == 63
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 63
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 63
               pshufd $0x55, rb3, rB1
            #elif KB == 63
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 63
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 63
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 63
               movaps 1888(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 63
               movaps 1904(pA), rA1
            #endif
         #endif
         #if KB > 63
            #if KB > 64
               movaps 896(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 64
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 64
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 64
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 64
               movaps 1920(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 64
               movaps 1936(pA), rA1
            #endif
         #endif
         #if KB > 64
            #if KB > 65
               movaps 912(pB), rb3
            #elif KB == 65
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 65
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 65
               pshufd $0x00, rb3, rB0
            #elif KB == 65
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 65
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 65
               pshufd $0x55, rb3, rB1
            #elif KB == 65
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 65
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 65
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 65
               movaps 1952(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 65
               movaps 1968(pA), rA1
            #endif
         #endif
         #if KB > 65
            #if KB > 66
               movaps 928(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 66
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 66
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 66
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 66
               movaps 1984(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 66
               movaps 2000(pA), rA1
            #endif
         #endif
         #if KB > 66
            #if KB > 67
               movaps 944(pB), rb3
            #elif KB == 67
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 67
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 67
               pshufd $0x00, rb3, rB0
            #elif KB == 67
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 67
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 67
               pshufd $0x55, rb3, rB1
            #elif KB == 67
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 67
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 67
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 67
               movaps 2016(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 67
               movaps 2032(pA), rA1
            #endif
         #endif
         #if KB > 67
            #if KB > 68
               movaps 960(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 68
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 68
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 68
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 68
               movaps 2048(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 68
               movaps 2064(pA), rA1
            #endif
         #endif
         #if KB > 68
            #if KB > 69
               movaps 976(pB), rb3
            #elif KB == 69
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 69
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 69
               pshufd $0x00, rb3, rB0
            #elif KB == 69
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 69
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 69
               pshufd $0x55, rb3, rB1
            #elif KB == 69
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 69
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 69
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 69
               movaps 2080(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 69
               movaps 2096(pA), rA1
            #endif
         #endif
         #if KB > 69
            #if KB > 70
               movaps 992(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 70
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 70
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 70
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 70
               movaps 2112(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 70
               movaps 2128(pA), rA1
            #endif
         #endif
         #if KB > 70
            #if KB > 71
               movaps 1008(pB), rb3
            #elif KB == 71
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 71
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 71
               pshufd $0x00, rb3, rB0
            #elif KB == 71
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 71
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 71
               pshufd $0x55, rb3, rB1
            #elif KB == 71
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 71
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 71
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 71
               movaps 2144(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 71
               movaps 2160(pA), rA1
            #endif
         #endif
         #if KB > 71
            #if KB > 72
               movaps 1024(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 72
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 72
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 72
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 72
               movaps 2176(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 72
               movaps 2192(pA), rA1
            #endif
         #endif
         #if KB > 72
            #if KB > 73
               movaps 1040(pB), rb3
            #elif KB == 73
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 73
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 73
               pshufd $0x00, rb3, rB0
            #elif KB == 73
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 73
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 73
               pshufd $0x55, rb3, rB1
            #elif KB == 73
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 73
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 73
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 73
               movaps 2208(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 73
               movaps 2224(pA), rA1
            #endif
         #endif
         #if KB > 73
            #if KB > 74
               movaps 1056(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 74
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 74
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 74
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 74
               movaps 2240(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 74
               movaps 2256(pA), rA1
            #endif
         #endif
         #if KB > 74
            #if KB > 75
               movaps 1072(pB), rb3
            #elif KB == 75
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 75
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 75
               pshufd $0x00, rb3, rB0
            #elif KB == 75
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 75
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 75
               pshufd $0x55, rb3, rB1
            #elif KB == 75
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 75
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 75
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 75
               movaps 2272(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 75
               movaps 2288(pA), rA1
            #endif
         #endif
         #if KB > 75
            #if KB > 76
               movaps 1088(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 76
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 76
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 76
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 76
               movaps 2304(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 76
               movaps 2320(pA), rA1
            #endif
         #endif
         #if KB > 76
            #if KB > 77
               movaps 1104(pB), rb3
            #elif KB == 77
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 77
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 77
               pshufd $0x00, rb3, rB0
            #elif KB == 77
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 77
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 77
               pshufd $0x55, rb3, rB1
            #elif KB == 77
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 77
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 77
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 77
               movaps 2336(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 77
               movaps 2352(pA), rA1
            #endif
         #endif
         #if KB > 77
            #if KB > 78
               movaps 1120(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 78
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 78
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 78
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 78
               movaps 2368(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 78
               movaps 2384(pA), rA1
            #endif
         #endif
         #if KB > 78
            #if KB > 79
               movaps 1136(pB), rb3
            #elif KB == 79
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 79
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 79
               pshufd $0x00, rb3, rB0
            #elif KB == 79
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 79
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 79
               pshufd $0x55, rb3, rB1
            #elif KB == 79
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 79
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 79
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 79
               movaps 2400(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 79
               movaps 2416(pA), rA1
            #endif
         #endif
         #if KB > 79
            #if KB > 80
               movaps 1152(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 80
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 80
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 80
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 80
               movaps 2432(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 80
               movaps 2448(pA), rA1
            #endif
         #endif
         #if KB > 80
            #if KB > 81
               movaps 1168(pB), rb3
            #elif KB == 81
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 81
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 81
               pshufd $0x00, rb3, rB0
            #elif KB == 81
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 81
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 81
               pshufd $0x55, rb3, rB1
            #elif KB == 81
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 81
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 81
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 81
               movaps 2464(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 81
               movaps 2480(pA), rA1
            #endif
         #endif
         #if KB > 81
            #if KB > 82
               movaps 1184(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 82
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 82
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 82
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 82
               movaps 2496(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 82
               movaps 2512(pA), rA1
            #endif
         #endif
         #if KB > 82
            #if KB > 83
               movaps 1200(pB), rb3
            #elif KB == 83
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 83
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 83
               pshufd $0x00, rb3, rB0
            #elif KB == 83
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 83
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 83
               pshufd $0x55, rb3, rB1
            #elif KB == 83
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 83
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 83
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 83
               movaps 2528(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 83
               movaps 2544(pA), rA1
            #endif
         #endif
         #if KB > 83
            #if KB > 84
               movaps 1216(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 84
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 84
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 84
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 84
               movaps 2560(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 84
               movaps 2576(pA), rA1
            #endif
         #endif
         #if KB > 84
            #if KB > 85
               movaps 1232(pB), rb3
            #elif KB == 85
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 85
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 85
               pshufd $0x00, rb3, rB0
            #elif KB == 85
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 85
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 85
               pshufd $0x55, rb3, rB1
            #elif KB == 85
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 85
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 85
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 85
               movaps 2592(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 85
               movaps 2608(pA), rA1
            #endif
         #endif
         #if KB > 85
            #if KB > 86
               movaps 1248(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 86
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 86
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 86
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 86
               movaps 2624(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 86
               movaps 2640(pA), rA1
            #endif
         #endif
         #if KB > 86
            #if KB > 87
               movaps 1264(pB), rb3
            #elif KB == 87
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 87
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 87
               pshufd $0x00, rb3, rB0
            #elif KB == 87
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 87
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 87
               pshufd $0x55, rb3, rB1
            #elif KB == 87
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 87
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 87
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 87
               movaps 2656(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 87
               movaps 2672(pA), rA1
            #endif
         #endif
         #if KB > 87
            #if KB > 88
               movaps 1280(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 88
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 88
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 88
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 88
               movaps 2688(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 88
               movaps 2704(pA), rA1
            #endif
         #endif
         #if KB > 88
            #if KB > 89
               movaps 1296(pB), rb3
            #elif KB == 89
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 89
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 89
               pshufd $0x00, rb3, rB0
            #elif KB == 89
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 89
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 89
               pshufd $0x55, rb3, rB1
            #elif KB == 89
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 89
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 89
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 89
               movaps 2720(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 89
               movaps 2736(pA), rA1
            #endif
         #endif
         #if KB > 89
            #if KB > 90
               movaps 1312(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 90
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 90
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 90
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 90
               movaps 2752(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 90
               movaps 2768(pA), rA1
            #endif
         #endif
         #if KB > 90
            #if KB > 91
               movaps 1328(pB), rb3
            #elif KB == 91
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 91
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 91
               pshufd $0x00, rb3, rB0
            #elif KB == 91
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 91
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 91
               pshufd $0x55, rb3, rB1
            #elif KB == 91
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 91
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 91
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 91
               movaps 2784(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 91
               movaps 2800(pA), rA1
            #endif
         #endif
         #if KB > 91
            #if KB > 92
               movaps 1344(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 92
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 92
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 92
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 92
               movaps 2816(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 92
               movaps 2832(pA), rA1
            #endif
         #endif
         #if KB > 92
            #if KB > 93
               movaps 1360(pB), rb3
            #elif KB == 93
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 93
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 93
               pshufd $0x00, rb3, rB0
            #elif KB == 93
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 93
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 93
               pshufd $0x55, rb3, rB1
            #elif KB == 93
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 93
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 93
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 93
               movaps 2848(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 93
               movaps 2864(pA), rA1
            #endif
         #endif
         #if KB > 93
            #if KB > 94
               movaps 1376(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 94
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 94
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 94
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 94
               movaps 2880(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 94
               movaps 2896(pA), rA1
            #endif
         #endif
         #if KB > 94
            #if KB > 95
               movaps 1392(pB), rb3
            #elif KB == 95
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 95
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 95
               pshufd $0x00, rb3, rB0
            #elif KB == 95
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 95
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 95
               pshufd $0x55, rb3, rB1
            #elif KB == 95
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 95
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 95
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 95
               movaps 2912(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 95
               movaps 2928(pA), rA1
            #endif
         #endif
         #if KB > 95
            #if KB > 96
               movaps 1408(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 96
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 96
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 96
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 96
               movaps 2944(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 96
               movaps 2960(pA), rA1
            #endif
         #endif
         #if KB > 96
            #if KB > 97
               movaps 1424(pB), rb3
            #elif KB == 97
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 97
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 97
               pshufd $0x00, rb3, rB0
            #elif KB == 97
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 97
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 97
               pshufd $0x55, rb3, rB1
            #elif KB == 97
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 97
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 97
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 97
               movaps 2976(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 97
               movaps 2992(pA), rA1
            #endif
         #endif
         #if KB > 97
            #if KB > 98
               movaps 1440(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 98
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 98
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 98
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 98
               movaps 3008(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 98
               movaps 3024(pA), rA1
            #endif
         #endif
         #if KB > 98
            #if KB > 99
               movaps 1456(pB), rb3
            #elif KB == 99
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 99
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 99
               pshufd $0x00, rb3, rB0
            #elif KB == 99
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 99
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 99
               pshufd $0x55, rb3, rB1
            #elif KB == 99
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 99
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 99
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 99
               movaps 3040(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 99
               movaps 3056(pA), rA1
            #endif
         #endif
         #if KB > 99
            #if KB > 100
               movaps 1472(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 100
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 100
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 100
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 100
               movaps 3072(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 100
               movaps 3088(pA), rA1
            #endif
         #endif
         #if KB > 100
            #if KB > 101
               movaps 1488(pB), rb3
            #elif KB == 101
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 101
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 101
               pshufd $0x00, rb3, rB0
            #elif KB == 101
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 101
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 101
               pshufd $0x55, rb3, rB1
            #elif KB == 101
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 101
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 101
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 101
               movaps 3104(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 101
               movaps 3120(pA), rA1
            #endif
         #endif
         #if KB > 101
            #if KB > 102
               movaps 1504(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 102
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 102
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 102
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 102
               movaps 3136(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 102
               movaps 3152(pA), rA1
            #endif
         #endif
         #if KB > 102
            #if KB > 103
               movaps 1520(pB), rb3
            #elif KB == 103
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 103
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 103
               pshufd $0x00, rb3, rB0
            #elif KB == 103
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 103
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 103
               pshufd $0x55, rb3, rB1
            #elif KB == 103
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 103
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 103
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 103
               movaps 3168(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 103
               movaps 3184(pA), rA1
            #endif
         #endif
         #if KB > 103
            #if KB > 104
               movaps 1536(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 104
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 104
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 104
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 104
               movaps 3200(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 104
               movaps 3216(pA), rA1
            #endif
         #endif
         #if KB > 104
            #if KB > 105
               movaps 1552(pB), rb3
            #elif KB == 105
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 105
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 105
               pshufd $0x00, rb3, rB0
            #elif KB == 105
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 105
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 105
               pshufd $0x55, rb3, rB1
            #elif KB == 105
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 105
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 105
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 105
               movaps 3232(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 105
               movaps 3248(pA), rA1
            #endif
         #endif
         #if KB > 105
            #if KB > 106
               movaps 1568(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 106
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 106
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 106
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 106
               movaps 3264(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 106
               movaps 3280(pA), rA1
            #endif
         #endif
         #if KB > 106
            #if KB > 107
               movaps 1584(pB), rb3
            #elif KB == 107
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 107
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 107
               pshufd $0x00, rb3, rB0
            #elif KB == 107
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 107
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 107
               pshufd $0x55, rb3, rB1
            #elif KB == 107
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 107
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 107
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 107
               movaps 3296(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 107
               movaps 3312(pA), rA1
            #endif
         #endif
         #if KB > 107
            #if KB > 108
               movaps 1600(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 108
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 108
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 108
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 108
               movaps 3328(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 108
               movaps 3344(pA), rA1
            #endif
         #endif
         #if KB > 108
            #if KB > 109
               movaps 1616(pB), rb3
            #elif KB == 109
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 109
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 109
               pshufd $0x00, rb3, rB0
            #elif KB == 109
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 109
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 109
               pshufd $0x55, rb3, rB1
            #elif KB == 109
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 109
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 109
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 109
               movaps 3360(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 109
               movaps 3376(pA), rA1
            #endif
         #endif
         #if KB > 109
            #if KB > 110
               movaps 1632(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 110
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 110
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 110
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 110
               movaps 3392(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 110
               movaps 3408(pA), rA1
            #endif
         #endif
         #if KB > 110
            #if KB > 111
               movaps 1648(pB), rb3
            #elif KB == 111
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 111
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 111
               pshufd $0x00, rb3, rB0
            #elif KB == 111
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 111
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 111
               pshufd $0x55, rb3, rB1
            #elif KB == 111
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 111
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 111
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 111
               movaps 3424(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 111
               movaps 3440(pA), rA1
            #endif
         #endif
         #if KB > 111
            #if KB > 112
               movaps 1664(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 112
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 112
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 112
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 112
               movaps 3456(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 112
               movaps 3472(pA), rA1
            #endif
         #endif
         #if KB > 112
            #if KB > 113
               movaps 1680(pB), rb3
            #elif KB == 113
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 113
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 113
               pshufd $0x00, rb3, rB0
            #elif KB == 113
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 113
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 113
               pshufd $0x55, rb3, rB1
            #elif KB == 113
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 113
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 113
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 113
               movaps 3488(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 113
               movaps 3504(pA), rA1
            #endif
         #endif
         #if KB > 113
            #if KB > 114
               movaps 1696(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 114
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 114
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 114
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 114
               movaps 3520(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 114
               movaps 3536(pA), rA1
            #endif
         #endif
         #if KB > 114
            #if KB > 115
               movaps 1712(pB), rb3
            #elif KB == 115
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 115
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 115
               pshufd $0x00, rb3, rB0
            #elif KB == 115
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 115
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 115
               pshufd $0x55, rb3, rB1
            #elif KB == 115
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 115
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 115
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 115
               movaps 3552(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 115
               movaps 3568(pA), rA1
            #endif
         #endif
         #if KB > 115
            #if KB > 116
               movaps 1728(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 116
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 116
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 116
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 116
               movaps 3584(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 116
               movaps 3600(pA), rA1
            #endif
         #endif
         #if KB > 116
            #if KB > 117
               movaps 1744(pB), rb3
            #elif KB == 117
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 117
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 117
               pshufd $0x00, rb3, rB0
            #elif KB == 117
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 117
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 117
               pshufd $0x55, rb3, rB1
            #elif KB == 117
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 117
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 117
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 117
               movaps 3616(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 117
               movaps 3632(pA), rA1
            #endif
         #endif
         #if KB > 117
            #if KB > 118
               movaps 1760(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 118
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 118
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 118
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 118
               movaps 3648(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 118
               movaps 3664(pA), rA1
            #endif
         #endif
         #if KB > 118
            #if KB > 119
               movaps 1776(pB), rb3
            #elif KB == 119
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 119
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 119
               pshufd $0x00, rb3, rB0
            #elif KB == 119
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 119
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 119
               pshufd $0x55, rb3, rB1
            #elif KB == 119
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 119
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 119
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 119
               movaps 3680(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 119
               movaps 3696(pA), rA1
            #endif
         #endif
         #if KB > 119
            #if KB > 120
               movaps 1792(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 120
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 120
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 120
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 120
               movaps 3712(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 120
               movaps 3728(pA), rA1
            #endif
         #endif
         #if KB > 120
            #if KB > 121
               movaps 1808(pB), rb3
            #elif KB == 121
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 121
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 121
               pshufd $0x00, rb3, rB0
            #elif KB == 121
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 121
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 121
               pshufd $0x55, rb3, rB1
            #elif KB == 121
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 121
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 121
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 121
               movaps 3744(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 121
               movaps 3760(pA), rA1
            #endif
         #endif
         #if KB > 121
            #if KB > 122
               movaps 1824(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 122
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 122
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 122
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 122
               movaps 3776(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 122
               movaps 3792(pA), rA1
            #endif
         #endif
         #if KB > 122
            #if KB > 123
               movaps 1840(pB), rb3
            #elif KB == 123
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 123
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 123
               pshufd $0x00, rb3, rB0
            #elif KB == 123
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 123
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 123
               pshufd $0x55, rb3, rB1
            #elif KB == 123
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 123
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 123
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 123
               movaps 3808(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 123
               movaps 3824(pA), rA1
            #endif
         #endif
         #if KB > 123
            #if KB > 124
               movaps 1856(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 124
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 124
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 124
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 124
               movaps 3840(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 124
               movaps 3856(pA), rA1
            #endif
         #endif
         #if KB > 124
            #if KB > 125
               movaps 1872(pB), rb3
            #elif KB == 125
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 125
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 125
               pshufd $0x00, rb3, rB0
            #elif KB == 125
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 125
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 125
               pshufd $0x55, rb3, rB1
            #elif KB == 125
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 125
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 125
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 125
               movaps 3872(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 125
               movaps 3888(pA), rA1
            #endif
         #endif
         #if KB > 125
            #if KB > 126
               movaps 1888(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 126
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 126
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 126
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 126
               movaps 3904(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 126
               movaps 3920(pA), rA1
            #endif
         #endif
         #if KB > 126
            #if KB > 127
               movaps 1904(pB), rb3
            #elif KB == 127
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 127
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 127
               pshufd $0x00, rb3, rB0
            #elif KB == 127
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 127
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 127
               pshufd $0x55, rb3, rB1
            #elif KB == 127
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 127
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 127
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 127
               movaps 3936(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 127
               movaps 3952(pA), rA1
            #endif
         #endif
         #if KB > 127
            #if KB > 128
               movaps 1920(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 128
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 128
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 128
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 128
               movaps 3968(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 128
               movaps 3984(pA), rA1
            #endif
         #endif
         #if KB > 128
            #if KB > 129
               movaps 1936(pB), rb3
            #elif KB == 129
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 129
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 129
               pshufd $0x00, rb3, rB0
            #elif KB == 129
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 129
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 129
               pshufd $0x55, rb3, rB1
            #elif KB == 129
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 129
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 129
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 129
               movaps 4000(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 129
               movaps 4016(pA), rA1
            #endif
         #endif
         #if KB > 129
            #if KB > 130
               movaps 1952(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 130
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 130
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 130
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 130
               movaps 4032(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 130
               movaps 4048(pA), rA1
            #endif
         #endif
         #if KB > 130
            #if KB > 131
               movaps 1968(pB), rb3
            #elif KB == 131
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 131
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 131
               pshufd $0x00, rb3, rB0
            #elif KB == 131
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 131
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 131
               pshufd $0x55, rb3, rB1
            #elif KB == 131
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 131
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 131
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 131
               movaps 4064(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 131
               movaps 4080(pA), rA1
            #endif
         #endif
         #if KB > 131
            #if KB > 132
               movaps 1984(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 132
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 132
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 132
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 132
               movaps 4096(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 132
               movaps 4112(pA), rA1
            #endif
         #endif
         #if KB > 132
            #if KB > 133
               movaps 2000(pB), rb3
            #elif KB == 133
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 133
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 133
               pshufd $0x00, rb3, rB0
            #elif KB == 133
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 133
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 133
               pshufd $0x55, rb3, rB1
            #elif KB == 133
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 133
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 133
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 133
               movaps 4128(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 133
               movaps 4144(pA), rA1
            #endif
         #endif
         #if KB > 133
            #if KB > 134
               movaps 2016(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 134
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 134
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 134
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 134
               movaps 4160(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 134
               movaps 4176(pA), rA1
            #endif
         #endif
         #if KB > 134
            #if KB > 135
               movaps 2032(pB), rb3
            #elif KB == 135
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 135
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 135
               pshufd $0x00, rb3, rB0
            #elif KB == 135
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 135
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 135
               pshufd $0x55, rb3, rB1
            #elif KB == 135
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 135
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 135
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 135
               movaps 4192(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 135
               movaps 4208(pA), rA1
            #endif
         #endif
         #if KB > 135
            #if KB > 136
               movaps 2048(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 136
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 136
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 136
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 136
               movaps 4224(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 136
               movaps 4240(pA), rA1
            #endif
         #endif
         #if KB > 136
            #if KB > 137
               movaps 2064(pB), rb3
            #elif KB == 137
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 137
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 137
               pshufd $0x00, rb3, rB0
            #elif KB == 137
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 137
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 137
               pshufd $0x55, rb3, rB1
            #elif KB == 137
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 137
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 137
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 137
               movaps 4256(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 137
               movaps 4272(pA), rA1
            #endif
         #endif
         #if KB > 137
            #if KB > 138
               movaps 2080(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 138
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 138
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 138
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 138
               movaps 4288(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 138
               movaps 4304(pA), rA1
            #endif
         #endif
         #if KB > 138
            #if KB > 139
               movaps 2096(pB), rb3
            #elif KB == 139
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 139
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 139
               pshufd $0x00, rb3, rB0
            #elif KB == 139
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 139
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 139
               pshufd $0x55, rb3, rB1
            #elif KB == 139
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 139
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 139
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 139
               movaps 4320(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 139
               movaps 4336(pA), rA1
            #endif
         #endif
         #if KB > 139
            #if KB > 140
               movaps 2112(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 140
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 140
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 140
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 140
               movaps 4352(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 140
               movaps 4368(pA), rA1
            #endif
         #endif
         #if KB > 140
            #if KB > 141
               movaps 2128(pB), rb3
            #elif KB == 141
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 141
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 141
               pshufd $0x00, rb3, rB0
            #elif KB == 141
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 141
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 141
               pshufd $0x55, rb3, rB1
            #elif KB == 141
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 141
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 141
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 141
               movaps 4384(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 141
               movaps 4400(pA), rA1
            #endif
         #endif
         #if KB > 141
            #if KB > 142
               movaps 2144(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 142
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 142
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 142
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 142
               movaps 4416(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 142
               movaps 4432(pA), rA1
            #endif
         #endif
         #if KB > 142
            #if KB > 143
               movaps 2160(pB), rb3
            #elif KB == 143
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 143
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 143
               pshufd $0x00, rb3, rB0
            #elif KB == 143
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 143
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 143
               pshufd $0x55, rb3, rB1
            #elif KB == 143
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 143
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 143
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 143
               movaps 4448(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 143
               movaps 4464(pA), rA1
            #endif
         #endif
         #if KB > 143
            #if KB > 144
               movaps 2176(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 144
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 144
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 144
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 144
               movaps 4480(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 144
               movaps 4496(pA), rA1
            #endif
         #endif
         #if KB > 144
            #if KB > 145
               movaps 2192(pB), rb3
            #elif KB == 145
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 145
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 145
               pshufd $0x00, rb3, rB0
            #elif KB == 145
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 145
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 145
               pshufd $0x55, rb3, rB1
            #elif KB == 145
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 145
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 145
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 145
               movaps 4512(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 145
               movaps 4528(pA), rA1
            #endif
         #endif
         #if KB > 145
            #if KB > 146
               movaps 2208(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 146
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 146
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 146
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 146
               movaps 4544(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 146
               movaps 4560(pA), rA1
            #endif
         #endif
         #if KB > 146
            #if KB > 147
               movaps 2224(pB), rb3
            #elif KB == 147
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 147
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 147
               pshufd $0x00, rb3, rB0
            #elif KB == 147
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 147
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 147
               pshufd $0x55, rb3, rB1
            #elif KB == 147
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 147
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 147
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 147
               movaps 4576(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 147
               movaps 4592(pA), rA1
            #endif
         #endif
         #if KB > 147
            #if KB > 148
               movaps 2240(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 148
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 148
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 148
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 148
               movaps 4608(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 148
               movaps 4624(pA), rA1
            #endif
         #endif
         #if KB > 148
            #if KB > 149
               movaps 2256(pB), rb3
            #elif KB == 149
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 149
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 149
               pshufd $0x00, rb3, rB0
            #elif KB == 149
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 149
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 149
               pshufd $0x55, rb3, rB1
            #elif KB == 149
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 149
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 149
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 149
               movaps 4640(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 149
               movaps 4656(pA), rA1
            #endif
         #endif
         #if KB > 149
            #if KB > 150
               movaps 2272(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 150
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 150
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 150
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 150
               movaps 4672(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 150
               movaps 4688(pA), rA1
            #endif
         #endif
         #if KB > 150
            #if KB > 151
               movaps 2288(pB), rb3
            #elif KB == 151
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 151
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 151
               pshufd $0x00, rb3, rB0
            #elif KB == 151
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 151
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 151
               pshufd $0x55, rb3, rB1
            #elif KB == 151
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 151
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 151
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 151
               movaps 4704(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 151
               movaps 4720(pA), rA1
            #endif
         #endif
         #if KB > 151
            #if KB > 152
               movaps 2304(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 152
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 152
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 152
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 152
               movaps 4736(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 152
               movaps 4752(pA), rA1
            #endif
         #endif
         #if KB > 152
            #if KB > 153
               movaps 2320(pB), rb3
            #elif KB == 153
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 153
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 153
               pshufd $0x00, rb3, rB0
            #elif KB == 153
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 153
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 153
               pshufd $0x55, rb3, rB1
            #elif KB == 153
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 153
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 153
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 153
               movaps 4768(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 153
               movaps 4784(pA), rA1
            #endif
         #endif
         #if KB > 153
            #if KB > 154
               movaps 2336(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 154
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 154
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 154
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 154
               movaps 4800(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 154
               movaps 4816(pA), rA1
            #endif
         #endif
         #if KB > 154
            #if KB > 155
               movaps 2352(pB), rb3
            #elif KB == 155
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 155
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 155
               pshufd $0x00, rb3, rB0
            #elif KB == 155
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 155
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 155
               pshufd $0x55, rb3, rB1
            #elif KB == 155
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 155
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 155
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 155
               movaps 4832(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 155
               movaps 4848(pA), rA1
            #endif
         #endif
         #if KB > 155
            #if KB > 156
               movaps 2368(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 156
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 156
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 156
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 156
               movaps 4864(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 156
               movaps 4880(pA), rA1
            #endif
         #endif
         #if KB > 156
            #if KB > 157
               movaps 2384(pB), rb3
            #elif KB == 157
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 157
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 157
               pshufd $0x00, rb3, rB0
            #elif KB == 157
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 157
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 157
               pshufd $0x55, rb3, rB1
            #elif KB == 157
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 157
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 157
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 157
               movaps 4896(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 157
               movaps 4912(pA), rA1
            #endif
         #endif
         #if KB > 157
            #if KB > 158
               movaps 2400(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 158
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 158
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 158
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 158
               movaps 4928(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 158
               movaps 4944(pA), rA1
            #endif
         #endif
         #if KB > 158
            #if KB > 159
               movaps 2416(pB), rb3
            #elif KB == 159
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 159
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 159
               pshufd $0x00, rb3, rB0
            #elif KB == 159
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 159
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 159
               pshufd $0x55, rb3, rB1
            #elif KB == 159
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 159
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 159
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 159
               movaps 4960(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 159
               movaps 4976(pA), rA1
            #endif
         #endif
         #if KB > 159
            #if KB > 160
               movaps 2432(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 160
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 160
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 160
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 160
               movaps 4992(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 160
               movaps 5008(pA), rA1
            #endif
         #endif
         #if KB > 160
            #if KB > 161
               movaps 2448(pB), rb3
            #elif KB == 161
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 161
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 161
               pshufd $0x00, rb3, rB0
            #elif KB == 161
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 161
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 161
               pshufd $0x55, rb3, rB1
            #elif KB == 161
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 161
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 161
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 161
               movaps 5024(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 161
               movaps 5040(pA), rA1
            #endif
         #endif
         #if KB > 161
            #if KB > 162
               movaps 2464(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 162
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 162
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 162
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 162
               movaps 5056(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 162
               movaps 5072(pA), rA1
            #endif
         #endif
         #if KB > 162
            #if KB > 163
               movaps 2480(pB), rb3
            #elif KB == 163
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 163
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 163
               pshufd $0x00, rb3, rB0
            #elif KB == 163
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 163
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 163
               pshufd $0x55, rb3, rB1
            #elif KB == 163
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 163
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 163
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 163
               movaps 5088(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 163
               movaps 5104(pA), rA1
            #endif
         #endif
         #if KB > 163
            #if KB > 164
               movaps 2496(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 164
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 164
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 164
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 164
               movaps 5120(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 164
               movaps 5136(pA), rA1
            #endif
         #endif
         #if KB > 164
            #if KB > 165
               movaps 2512(pB), rb3
            #elif KB == 165
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 165
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 165
               pshufd $0x00, rb3, rB0
            #elif KB == 165
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 165
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 165
               pshufd $0x55, rb3, rB1
            #elif KB == 165
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 165
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 165
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 165
               movaps 5152(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 165
               movaps 5168(pA), rA1
            #endif
         #endif
         #if KB > 165
            #if KB > 166
               movaps 2528(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 166
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 166
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 166
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 166
               movaps 5184(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 166
               movaps 5200(pA), rA1
            #endif
         #endif
         #if KB > 166
            #if KB > 167
               movaps 2544(pB), rb3
            #elif KB == 167
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 167
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 167
               pshufd $0x00, rb3, rB0
            #elif KB == 167
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 167
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 167
               pshufd $0x55, rb3, rB1
            #elif KB == 167
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 167
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 167
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 167
               movaps 5216(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 167
               movaps 5232(pA), rA1
            #endif
         #endif
         #if KB > 167
            #if KB > 168
               movaps 2560(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 168
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 168
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 168
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 168
               movaps 5248(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 168
               movaps 5264(pA), rA1
            #endif
         #endif
         #if KB > 168
            #if KB > 169
               movaps 2576(pB), rb3
            #elif KB == 169
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 169
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 169
               pshufd $0x00, rb3, rB0
            #elif KB == 169
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 169
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 169
               pshufd $0x55, rb3, rB1
            #elif KB == 169
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 169
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 169
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 169
               movaps 5280(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 169
               movaps 5296(pA), rA1
            #endif
         #endif
         #if KB > 169
            #if KB > 170
               movaps 2592(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 170
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 170
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 170
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 170
               movaps 5312(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 170
               movaps 5328(pA), rA1
            #endif
         #endif
         #if KB > 170
            #if KB > 171
               movaps 2608(pB), rb3
            #elif KB == 171
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 171
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 171
               pshufd $0x00, rb3, rB0
            #elif KB == 171
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 171
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 171
               pshufd $0x55, rb3, rB1
            #elif KB == 171
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 171
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 171
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 171
               movaps 5344(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 171
               movaps 5360(pA), rA1
            #endif
         #endif
         #if KB > 171
            #if KB > 172
               movaps 2624(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 172
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 172
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 172
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 172
               movaps 5376(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 172
               movaps 5392(pA), rA1
            #endif
         #endif
         #if KB > 172
            #if KB > 173
               movaps 2640(pB), rb3
            #elif KB == 173
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 173
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 173
               pshufd $0x00, rb3, rB0
            #elif KB == 173
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 173
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 173
               pshufd $0x55, rb3, rB1
            #elif KB == 173
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 173
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 173
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 173
               movaps 5408(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 173
               movaps 5424(pA), rA1
            #endif
         #endif
         #if KB > 173
            #if KB > 174
               movaps 2656(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 174
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 174
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 174
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 174
               movaps 5440(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 174
               movaps 5456(pA), rA1
            #endif
         #endif
         #if KB > 174
            #if KB > 175
               movaps 2672(pB), rb3
            #elif KB == 175
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 175
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 175
               pshufd $0x00, rb3, rB0
            #elif KB == 175
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 175
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 175
               pshufd $0x55, rb3, rB1
            #elif KB == 175
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 175
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 175
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 175
               movaps 5472(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 175
               movaps 5488(pA), rA1
            #endif
         #endif
         #if KB > 175
            #if KB > 176
               movaps 2688(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 176
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 176
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 176
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 176
               movaps 5504(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 176
               movaps 5520(pA), rA1
            #endif
         #endif
         #if KB > 176
            #if KB > 177
               movaps 2704(pB), rb3
            #elif KB == 177
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 177
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 177
               pshufd $0x00, rb3, rB0
            #elif KB == 177
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 177
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 177
               pshufd $0x55, rb3, rB1
            #elif KB == 177
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 177
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 177
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 177
               movaps 5536(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 177
               movaps 5552(pA), rA1
            #endif
         #endif
         #if KB > 177
            #if KB > 178
               movaps 2720(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 178
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 178
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 178
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 178
               movaps 5568(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 178
               movaps 5584(pA), rA1
            #endif
         #endif
         #if KB > 178
            #if KB > 179
               movaps 2736(pB), rb3
            #elif KB == 179
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 179
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 179
               pshufd $0x00, rb3, rB0
            #elif KB == 179
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 179
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 179
               pshufd $0x55, rb3, rB1
            #elif KB == 179
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 179
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 179
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 179
               movaps 5600(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 179
               movaps 5616(pA), rA1
            #endif
         #endif
         #if KB > 179
            #if KB > 180
               movaps 2752(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 180
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 180
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 180
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 180
               movaps 5632(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 180
               movaps 5648(pA), rA1
            #endif
         #endif
         #if KB > 180
            #if KB > 181
               movaps 2768(pB), rb3
            #elif KB == 181
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 181
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 181
               pshufd $0x00, rb3, rB0
            #elif KB == 181
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 181
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 181
               pshufd $0x55, rb3, rB1
            #elif KB == 181
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 181
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 181
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 181
               movaps 5664(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 181
               movaps 5680(pA), rA1
            #endif
         #endif
         #if KB > 181
            #if KB > 182
               movaps 2784(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 182
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 182
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 182
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 182
               movaps 5696(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 182
               movaps 5712(pA), rA1
            #endif
         #endif
         #if KB > 182
            #if KB > 183
               movaps 2800(pB), rb3
            #elif KB == 183
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 183
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 183
               pshufd $0x00, rb3, rB0
            #elif KB == 183
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 183
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 183
               pshufd $0x55, rb3, rB1
            #elif KB == 183
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 183
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 183
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 183
               movaps 5728(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 183
               movaps 5744(pA), rA1
            #endif
         #endif
         #if KB > 183
            #if KB > 184
               movaps 2816(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 184
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 184
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 184
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 184
               movaps 5760(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 184
               movaps 5776(pA), rA1
            #endif
         #endif
         #if KB > 184
            #if KB > 185
               movaps 2832(pB), rb3
            #elif KB == 185
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 185
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 185
               pshufd $0x00, rb3, rB0
            #elif KB == 185
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 185
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 185
               pshufd $0x55, rb3, rB1
            #elif KB == 185
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 185
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 185
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 185
               movaps 5792(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 185
               movaps 5808(pA), rA1
            #endif
         #endif
         #if KB > 185
            #if KB > 186
               movaps 2848(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 186
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 186
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 186
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 186
               movaps 5824(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 186
               movaps 5840(pA), rA1
            #endif
         #endif
         #if KB > 186
            #if KB > 187
               movaps 2864(pB), rb3
            #elif KB == 187
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 187
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 187
               pshufd $0x00, rb3, rB0
            #elif KB == 187
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 187
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 187
               pshufd $0x55, rb3, rB1
            #elif KB == 187
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 187
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 187
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 187
               movaps 5856(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 187
               movaps 5872(pA), rA1
            #endif
         #endif
         #if KB > 187
            #if KB > 188
               movaps 2880(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 188
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 188
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 188
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 188
               movaps 5888(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 188
               movaps 5904(pA), rA1
            #endif
         #endif
         #if KB > 188
            #if KB > 189
               movaps 2896(pB), rb3
            #elif KB == 189
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 189
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 189
               pshufd $0x00, rb3, rB0
            #elif KB == 189
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 189
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 189
               pshufd $0x55, rb3, rB1
            #elif KB == 189
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 189
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 189
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 189
               movaps 5920(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 189
               movaps 5936(pA), rA1
            #endif
         #endif
         #if KB > 189
            #if KB > 190
               movaps 2912(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 190
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 190
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 190
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 190
               movaps 5952(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 190
               movaps 5968(pA), rA1
            #endif
         #endif
         #if KB > 190
            #if KB > 191
               movaps 2928(pB), rb3
            #elif KB == 191
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 191
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 191
               pshufd $0x00, rb3, rB0
            #elif KB == 191
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 191
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 191
               pshufd $0x55, rb3, rB1
            #elif KB == 191
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 191
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 191
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 191
               movaps 5984(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 191
               movaps 6000(pA), rA1
            #endif
         #endif
         #if KB > 191
            #if KB > 192
               movaps 2944(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 192
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 192
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 192
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 192
               movaps 6016(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 192
               movaps 6032(pA), rA1
            #endif
         #endif
         #if KB > 192
            #if KB > 193
               movaps 2960(pB), rb3
            #elif KB == 193
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 193
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 193
               pshufd $0x00, rb3, rB0
            #elif KB == 193
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 193
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 193
               pshufd $0x55, rb3, rB1
            #elif KB == 193
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 193
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 193
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 193
               movaps 6048(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 193
               movaps 6064(pA), rA1
            #endif
         #endif
         #if KB > 193
            #if KB > 194
               movaps 2976(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 194
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 194
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 194
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 194
               movaps 6080(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 194
               movaps 6096(pA), rA1
            #endif
         #endif
         #if KB > 194
            #if KB > 195
               movaps 2992(pB), rb3
            #elif KB == 195
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 195
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 195
               pshufd $0x00, rb3, rB0
            #elif KB == 195
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 195
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 195
               pshufd $0x55, rb3, rB1
            #elif KB == 195
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 195
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 195
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 195
               movaps 6112(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 195
               movaps 6128(pA), rA1
            #endif
         #endif
         #if KB > 195
            #if KB > 196
               movaps 3008(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 196
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 196
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 196
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 196
               movaps 6144(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 196
               movaps 6160(pA), rA1
            #endif
         #endif
         #if KB > 196
            #if KB > 197
               movaps 3024(pB), rb3
            #elif KB == 197
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 197
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 197
               pshufd $0x00, rb3, rB0
            #elif KB == 197
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 197
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 197
               pshufd $0x55, rb3, rB1
            #elif KB == 197
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 197
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 197
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 197
               movaps 6176(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 197
               movaps 6192(pA), rA1
            #endif
         #endif
         #if KB > 197
            #if KB > 198
               movaps 3040(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 198
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 198
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 198
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 198
               movaps 6208(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 198
               movaps 6224(pA), rA1
            #endif
         #endif
         #if KB > 198
            #if KB > 199
               movaps 3056(pB), rb3
            #elif KB == 199
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 199
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 199
               pshufd $0x00, rb3, rB0
            #elif KB == 199
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 199
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 199
               pshufd $0x55, rb3, rB1
            #elif KB == 199
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 199
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 199
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 199
               movaps 6240(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 199
               movaps 6256(pA), rA1
            #endif
         #endif
         #if KB > 199
            #if KB > 200
               movaps 3072(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 200
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 200
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 200
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 200
               movaps 6272(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 200
               movaps 6288(pA), rA1
            #endif
         #endif
         #if KB > 200
            #if KB > 201
               movaps 3088(pB), rb3
            #elif KB == 201
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 201
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 201
               pshufd $0x00, rb3, rB0
            #elif KB == 201
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 201
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 201
               pshufd $0x55, rb3, rB1
            #elif KB == 201
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 201
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 201
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 201
               movaps 6304(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 201
               movaps 6320(pA), rA1
            #endif
         #endif
         #if KB > 201
            #if KB > 202
               movaps 3104(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 202
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 202
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 202
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 202
               movaps 6336(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 202
               movaps 6352(pA), rA1
            #endif
         #endif
         #if KB > 202
            #if KB > 203
               movaps 3120(pB), rb3
            #elif KB == 203
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 203
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 203
               pshufd $0x00, rb3, rB0
            #elif KB == 203
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 203
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 203
               pshufd $0x55, rb3, rB1
            #elif KB == 203
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 203
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 203
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 203
               movaps 6368(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 203
               movaps 6384(pA), rA1
            #endif
         #endif
         #if KB > 203
            #if KB > 204
               movaps 3136(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 204
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 204
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 204
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 204
               movaps 6400(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 204
               movaps 6416(pA), rA1
            #endif
         #endif
         #if KB > 204
            #if KB > 205
               movaps 3152(pB), rb3
            #elif KB == 205
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 205
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 205
               pshufd $0x00, rb3, rB0
            #elif KB == 205
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 205
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 205
               pshufd $0x55, rb3, rB1
            #elif KB == 205
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 205
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 205
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 205
               movaps 6432(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 205
               movaps 6448(pA), rA1
            #endif
         #endif
         #if KB > 205
            #if KB > 206
               movaps 3168(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 206
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 206
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 206
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 206
               movaps 6464(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 206
               movaps 6480(pA), rA1
            #endif
         #endif
         #if KB > 206
            #if KB > 207
               movaps 3184(pB), rb3
            #elif KB == 207
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 207
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 207
               pshufd $0x00, rb3, rB0
            #elif KB == 207
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 207
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 207
               pshufd $0x55, rb3, rB1
            #elif KB == 207
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 207
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 207
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 207
               movaps 6496(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 207
               movaps 6512(pA), rA1
            #endif
         #endif
         #if KB > 207
            #if KB > 208
               movaps 3200(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 208
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 208
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 208
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 208
               movaps 6528(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 208
               movaps 6544(pA), rA1
            #endif
         #endif
         #if KB > 208
            #if KB > 209
               movaps 3216(pB), rb3
            #elif KB == 209
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 209
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 209
               pshufd $0x00, rb3, rB0
            #elif KB == 209
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 209
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 209
               pshufd $0x55, rb3, rB1
            #elif KB == 209
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 209
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 209
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 209
               movaps 6560(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 209
               movaps 6576(pA), rA1
            #endif
         #endif
         #if KB > 209
            #if KB > 210
               movaps 3232(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 210
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 210
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 210
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 210
               movaps 6592(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 210
               movaps 6608(pA), rA1
            #endif
         #endif
         #if KB > 210
            #if KB > 211
               movaps 3248(pB), rb3
            #elif KB == 211
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 211
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 211
               pshufd $0x00, rb3, rB0
            #elif KB == 211
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 211
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 211
               pshufd $0x55, rb3, rB1
            #elif KB == 211
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 211
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 211
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 211
               movaps 6624(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 211
               movaps 6640(pA), rA1
            #endif
         #endif
         #if KB > 211
            #if KB > 212
               movaps 3264(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 212
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 212
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 212
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 212
               movaps 6656(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 212
               movaps 6672(pA), rA1
            #endif
         #endif
         #if KB > 212
            #if KB > 213
               movaps 3280(pB), rb3
            #elif KB == 213
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 213
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 213
               pshufd $0x00, rb3, rB0
            #elif KB == 213
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 213
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 213
               pshufd $0x55, rb3, rB1
            #elif KB == 213
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 213
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 213
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 213
               movaps 6688(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 213
               movaps 6704(pA), rA1
            #endif
         #endif
         #if KB > 213
            #if KB > 214
               movaps 3296(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 214
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 214
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 214
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 214
               movaps 6720(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 214
               movaps 6736(pA), rA1
            #endif
         #endif
         #if KB > 214
            #if KB > 215
               movaps 3312(pB), rb3
            #elif KB == 215
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 215
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 215
               pshufd $0x00, rb3, rB0
            #elif KB == 215
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 215
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 215
               pshufd $0x55, rb3, rB1
            #elif KB == 215
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 215
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 215
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 215
               movaps 6752(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 215
               movaps 6768(pA), rA1
            #endif
         #endif
         #if KB > 215
            #if KB > 216
               movaps 3328(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 216
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 216
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 216
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 216
               movaps 6784(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 216
               movaps 6800(pA), rA1
            #endif
         #endif
         #if KB > 216
            #if KB > 217
               movaps 3344(pB), rb3
            #elif KB == 217
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 217
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 217
               pshufd $0x00, rb3, rB0
            #elif KB == 217
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 217
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 217
               pshufd $0x55, rb3, rB1
            #elif KB == 217
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 217
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 217
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 217
               movaps 6816(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 217
               movaps 6832(pA), rA1
            #endif
         #endif
         #if KB > 217
            #if KB > 218
               movaps 3360(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 218
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 218
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 218
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 218
               movaps 6848(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 218
               movaps 6864(pA), rA1
            #endif
         #endif
         #if KB > 218
            #if KB > 219
               movaps 3376(pB), rb3
            #elif KB == 219
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 219
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 219
               pshufd $0x00, rb3, rB0
            #elif KB == 219
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 219
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 219
               pshufd $0x55, rb3, rB1
            #elif KB == 219
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 219
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 219
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 219
               movaps 6880(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 219
               movaps 6896(pA), rA1
            #endif
         #endif
         #if KB > 219
            #if KB > 220
               movaps 3392(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 220
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 220
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 220
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 220
               movaps 6912(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 220
               movaps 6928(pA), rA1
            #endif
         #endif
         #if KB > 220
            #if KB > 221
               movaps 3408(pB), rb3
            #elif KB == 221
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 221
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 221
               pshufd $0x00, rb3, rB0
            #elif KB == 221
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 221
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 221
               pshufd $0x55, rb3, rB1
            #elif KB == 221
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 221
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 221
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 221
               movaps 6944(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 221
               movaps 6960(pA), rA1
            #endif
         #endif
         #if KB > 221
            #if KB > 222
               movaps 3424(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 222
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 222
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 222
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 222
               movaps 6976(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 222
               movaps 6992(pA), rA1
            #endif
         #endif
         #if KB > 222
            #if KB > 223
               movaps 3440(pB), rb3
            #elif KB == 223
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 223
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 223
               pshufd $0x00, rb3, rB0
            #elif KB == 223
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 223
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 223
               pshufd $0x55, rb3, rB1
            #elif KB == 223
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 223
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 223
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 223
               movaps 7008(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 223
               movaps 7024(pA), rA1
            #endif
         #endif
         #if KB > 223
            #if KB > 224
               movaps 3456(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 224
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 224
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 224
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 224
               movaps 7040(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 224
               movaps 7056(pA), rA1
            #endif
         #endif
         #if KB > 224
            #if KB > 225
               movaps 3472(pB), rb3
            #elif KB == 225
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 225
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 225
               pshufd $0x00, rb3, rB0
            #elif KB == 225
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 225
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 225
               pshufd $0x55, rb3, rB1
            #elif KB == 225
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 225
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 225
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 225
               movaps 7072(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 225
               movaps 7088(pA), rA1
            #endif
         #endif
         #if KB > 225
            #if KB > 226
               movaps 3488(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 226
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 226
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 226
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 226
               movaps 7104(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 226
               movaps 7120(pA), rA1
            #endif
         #endif
         #if KB > 226
            #if KB > 227
               movaps 3504(pB), rb3
            #elif KB == 227
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 227
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 227
               pshufd $0x00, rb3, rB0
            #elif KB == 227
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 227
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 227
               pshufd $0x55, rb3, rB1
            #elif KB == 227
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 227
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 227
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 227
               movaps 7136(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 227
               movaps 7152(pA), rA1
            #endif
         #endif
         #if KB > 227
            #if KB > 228
               movaps 3520(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 228
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 228
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 228
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 228
               movaps 7168(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 228
               movaps 7184(pA), rA1
            #endif
         #endif
         #if KB > 228
            #if KB > 229
               movaps 3536(pB), rb3
            #elif KB == 229
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 229
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 229
               pshufd $0x00, rb3, rB0
            #elif KB == 229
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 229
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 229
               pshufd $0x55, rb3, rB1
            #elif KB == 229
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 229
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 229
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 229
               movaps 7200(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 229
               movaps 7216(pA), rA1
            #endif
         #endif
         #if KB > 229
            #if KB > 230
               movaps 3552(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 230
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 230
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 230
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 230
               movaps 7232(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 230
               movaps 7248(pA), rA1
            #endif
         #endif
         #if KB > 230
            #if KB > 231
               movaps 3568(pB), rb3
            #elif KB == 231
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 231
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 231
               pshufd $0x00, rb3, rB0
            #elif KB == 231
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 231
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 231
               pshufd $0x55, rb3, rB1
            #elif KB == 231
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 231
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 231
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 231
               movaps 7264(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 231
               movaps 7280(pA), rA1
            #endif
         #endif
         #if KB > 231
            #if KB > 232
               movaps 3584(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 232
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 232
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 232
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 232
               movaps 7296(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 232
               movaps 7312(pA), rA1
            #endif
         #endif
         #if KB > 232
            #if KB > 233
               movaps 3600(pB), rb3
            #elif KB == 233
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 233
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 233
               pshufd $0x00, rb3, rB0
            #elif KB == 233
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 233
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 233
               pshufd $0x55, rb3, rB1
            #elif KB == 233
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 233
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 233
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 233
               movaps 7328(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 233
               movaps 7344(pA), rA1
            #endif
         #endif
         #if KB > 233
            #if KB > 234
               movaps 3616(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 234
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 234
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 234
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 234
               movaps 7360(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 234
               movaps 7376(pA), rA1
            #endif
         #endif
         #if KB > 234
            #if KB > 235
               movaps 3632(pB), rb3
            #elif KB == 235
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 235
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 235
               pshufd $0x00, rb3, rB0
            #elif KB == 235
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 235
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 235
               pshufd $0x55, rb3, rB1
            #elif KB == 235
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 235
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 235
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 235
               movaps 7392(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 235
               movaps 7408(pA), rA1
            #endif
         #endif
         #if KB > 235
            #if KB > 236
               movaps 3648(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 236
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 236
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 236
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 236
               movaps 7424(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 236
               movaps 7440(pA), rA1
            #endif
         #endif
         #if KB > 236
            #if KB > 237
               movaps 3664(pB), rb3
            #elif KB == 237
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 237
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 237
               pshufd $0x00, rb3, rB0
            #elif KB == 237
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 237
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 237
               pshufd $0x55, rb3, rB1
            #elif KB == 237
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 237
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 237
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 237
               movaps 7456(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 237
               movaps 7472(pA), rA1
            #endif
         #endif
         #if KB > 237
            #if KB > 238
               movaps 3680(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 238
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 238
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 238
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 238
               movaps 7488(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 238
               movaps 7504(pA), rA1
            #endif
         #endif
         #if KB > 238
            #if KB > 239
               movaps 3696(pB), rb3
            #elif KB == 239
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 239
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 239
               pshufd $0x00, rb3, rB0
            #elif KB == 239
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 239
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 239
               pshufd $0x55, rb3, rB1
            #elif KB == 239
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 239
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 239
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 239
               movaps 7520(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 239
               movaps 7536(pA), rA1
            #endif
         #endif
         #if KB > 239
            #if KB > 240
               movaps 3712(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 240
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 240
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 240
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 240
               movaps 7552(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 240
               movaps 7568(pA), rA1
            #endif
         #endif
         #if KB > 240
            #if KB > 241
               movaps 3728(pB), rb3
            #elif KB == 241
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 241
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 241
               pshufd $0x00, rb3, rB0
            #elif KB == 241
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 241
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 241
               pshufd $0x55, rb3, rB1
            #elif KB == 241
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 241
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 241
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 241
               movaps 7584(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 241
               movaps 7600(pA), rA1
            #endif
         #endif
         #if KB > 241
            #if KB > 242
               movaps 3744(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 242
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 242
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 242
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 242
               movaps 7616(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 242
               movaps 7632(pA), rA1
            #endif
         #endif
         #if KB > 242
            #if KB > 243
               movaps 3760(pB), rb3
            #elif KB == 243
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 243
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 243
               pshufd $0x00, rb3, rB0
            #elif KB == 243
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 243
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 243
               pshufd $0x55, rb3, rB1
            #elif KB == 243
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 243
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 243
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 243
               movaps 7648(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 243
               movaps 7664(pA), rA1
            #endif
         #endif
         #if KB > 243
            #if KB > 244
               movaps 3776(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 244
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 244
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 244
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 244
               movaps 7680(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 244
               movaps 7696(pA), rA1
            #endif
         #endif
         #if KB > 244
            #if KB > 245
               movaps 3792(pB), rb3
            #elif KB == 245
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 245
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 245
               pshufd $0x00, rb3, rB0
            #elif KB == 245
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 245
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 245
               pshufd $0x55, rb3, rB1
            #elif KB == 245
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 245
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 245
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 245
               movaps 7712(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 245
               movaps 7728(pA), rA1
            #endif
         #endif
         #if KB > 245
            #if KB > 246
               movaps 3808(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 246
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 246
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 246
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 246
               movaps 7744(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 246
               movaps 7760(pA), rA1
            #endif
         #endif
         #if KB > 246
            #if KB > 247
               movaps 3824(pB), rb3
            #elif KB == 247
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 247
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 247
               pshufd $0x00, rb3, rB0
            #elif KB == 247
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 247
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 247
               pshufd $0x55, rb3, rB1
            #elif KB == 247
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 247
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 247
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 247
               movaps 7776(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 247
               movaps 7792(pA), rA1
            #endif
         #endif
         #if KB > 247
            #if KB > 248
               movaps 3840(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 248
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 248
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 248
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 248
               movaps 7808(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 248
               movaps 7824(pA), rA1
            #endif
         #endif
         #if KB > 248
            #if KB > 249
               movaps 3856(pB), rb3
            #elif KB == 249
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 249
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 249
               pshufd $0x00, rb3, rB0
            #elif KB == 249
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 249
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 249
               pshufd $0x55, rb3, rB1
            #elif KB == 249
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 249
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 249
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 249
               movaps 7840(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 249
               movaps 7856(pA), rA1
            #endif
         #endif
         #if KB > 249
            #if KB > 250
               movaps 3872(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 250
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 250
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 250
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 250
               movaps 7872(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 250
               movaps 7888(pA), rA1
            #endif
         #endif
         #if KB > 250
            #if KB > 251
               movaps 3888(pB), rb3
            #elif KB == 251
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 251
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 251
               pshufd $0x00, rb3, rB0
            #elif KB == 251
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 251
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 251
               pshufd $0x55, rb3, rB1
            #elif KB == 251
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 251
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 251
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 251
               movaps 7904(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 251
               movaps 7920(pA), rA1
            #endif
         #endif
         #if KB > 251
            #if KB > 252
               movaps 3904(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 252
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 252
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 252
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 252
               movaps 7936(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 252
               movaps 7952(pA), rA1
            #endif
         #endif
         #if KB > 252
            #if KB > 253
               movaps 3920(pB), rb3
            #elif KB == 253
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 253
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 253
               pshufd $0x00, rb3, rB0
            #elif KB == 253
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 253
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 253
               pshufd $0x55, rb3, rB1
            #elif KB == 253
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 253
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 253
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 253
               movaps 7968(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 253
               movaps 7984(pA), rA1
            #endif
         #endif
         #if KB > 253
            #if KB > 254
               movaps 3936(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 254
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 254
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 254
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 254
               movaps 8000(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 254
               movaps 8016(pA), rA1
            #endif
         #endif
         #if KB > 254
            #if KB > 255
               movaps 3952(pB), rb3
            #elif KB == 255
               pref (pf)
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rB3, rB3
            #if KB == 255
               add incPF, pf
            #endif
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 255
               pshufd $0x00, rb3, rB0
            #elif KB == 255
               pref (pfC)
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            #if KB == 255
               add incPF, pfC
            #endif
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 255
               pshufd $0x55, rb3, rB1
            #elif KB == 255
               pref (pfB)
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            #if KB == 255
               add incPF, pfB
            #endif
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 255
               pshufd $0xAA, rb3, rB2
            #endif

            mulps rB3, rA0
            addps rA0, rC03
            #if KB > 255
               movaps 8032(pA), rA0
            #endif
            mulps rA1, rB3
            addps rB3, rC13
            #if KB > 255
               movaps 8048(pA), rA1
            #endif
         #endif
         #if KB > 255
            #if KB > 256
               movaps 3968(pB), rB3
            #endif
            movaps rA0, rm0
            mulps rB0, rm0
            addps rm0, rC00
            shufps $0xFF, rb3, rb3
            movaps rA1, rm0
            mulps rB0, rm0
            addps rm0, rC10
            #if KB > 256
               pshufd $0x00, rB3, rB0
            #endif
            movaps rA0, rm0
            mulps rB1, rm0
            addps rm0, rC01
            movaps rA1, rm0
            mulps rB1, rm0
            addps rm0, rC11
            #if KB > 256
               pshufd $0x55, rB3, rB1
            #endif
            movaps rA0, rm0
            mulps rB2, rm0
            addps rm0, rC02
            movaps rA1, rm0
            mulps rB2, rm0
            addps rm0, rC12
            #if KB > 256
               pshufd $0xAA, rB3, rB2
            #endif

            mulps rb3, rA0
            addps rA0, rC03
            #if KB > 256
               movaps 8064(pA), rA0
            #endif
            mulps rA1, rb3
            addps rb3, rC13
            #if KB > 256
               movaps 8080(pA), rA1
            #endif
         #endif
         #if defined(BETA1) || defined(BETAN1)
            BETCOP (pC), rC00
            movaps rC00, (pC)
            BETCOP 16(pC), rC10
            movaps rC10, 16(pC)
            BETCOP 32(pC), rC01
            movaps rC01, 32(pC)
            BETCOP 48(pC), rC11
            movaps rC11, 48(pC)
            BETCOP 64(pC), rC02
            movaps rC02, 64(pC)
            BETCOP 80(pC), rC12
            movaps rC12, 80(pC)
            BETCOP 96(pC), rC03
            movaps rC03, 96(pC)
            BETCOP 112(pC), rC13
            movaps rC13, 112(pC)
         #else
            movaps rC00, (pC)
            movaps rC10, 16(pC)
            movaps rC01, 32(pC)
            movaps rC11, 48(pC)
            movaps rC02, 64(pC)
            movaps rC12, 80(pC)
            movaps rC03, 96(pC)
            movaps rC13, 112(pC)
         #endif
         sub $-128, pC
         add $KB*4*4, pB
      sub $1, nnu
      jnz MNLOOP
      mov nnu0, nnu
      mov pB0, pB
      add incAm, pA
   sub $1, nmu
   jnz MNLOOP

/* DONE: */
   movq    (%rsp), %rbp
   movq    8(%rsp), %rbx
   movq    16(%rsp), %r12
   movq    24(%rsp), %r13
   add $FSIZE, %rsp
   ret
