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
#define movapd movaps
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

#define rA0     %xmm0
#define rA1     %xmm1
#define rB0     %xmm2
#define rB1     %xmm3
#define rB2     %xmm4
#define rB3     %xmm5
#define rC00    %xmm6
#define rC10    %xmm7
#define rC01    %xmm8
#define rC11    %xmm9
#define rC02    %xmm10
#define rC12    %xmm11
#define rC03    %xmm12
#define rC13    %xmm13
#define rm0     %xmm14
#define unpckhpd movhlps
/* #define movddup pshufd $0x44, */
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
   #define BETCOP subpd
#else
   #define BETCOP addpd
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
   mov $8*4*4, incPF
/*
 * Extend range of small operands by starting at -128
 */
   sub $-128, pA
   sub $-128, pB
   mov $KB*4*8, incAm           /* incAm = KB*MU*size */
   movq pB, pB0

   ALIGN8
   .local MNLOOP
   MNLOOP:
/*
 *       Peel first iteration of K-loop to handle init of C to 0
 */
         movapd -128(pB), rB1
         movddup rB1, rB0
         movapd -128(pA), rC00

         movapd rC00, rC01
         mulpd rB0, rC00
         unpckhpd rB1, rB1
         movapd -112(pA), rC10
         movapd rC10, rC11
         mulpd rB0, rC10

         movapd -112(pB), rB3
         movddup rB3, rB2
         movapd rC01, rC02
         mulpd rB1, rC01
         unpckhpd rB3, rB3
         movapd rC11, rC12
         mulpd rB1, rC11
         #if KB > 1
            movapd -96(pB), rB1
         #else
            pref (pf)
         #endif

         #if KB > 1
            movapd -96(pA), rA0
         #else
            pref 64(pf)
         #endif
         movapd rC02, rC03
         mulpd rB2, rC02
         movapd rC12, rC13
         mulpd rB2, rC12
         #if KB > 1
            movddup rB1, rB0
         #else
            add incPF, pf
         #endif

         mulpd rB3, rC03
         prefC (pC)
         mulpd rB3, rC13
         prefC 64(pC)
/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 1
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd -80(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd -80(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 2
               movapd -64(pB), rB1
            #elif KB == 2
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 2
               movapd -64(pA), rA0
            #elif KB == 2
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 2
               movddup rB1, rB0
            #elif KB == 2
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 2
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd -48(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd -48(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 3
               movapd -32(pB), rB1
            #elif KB == 3
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 3
               movapd -32(pA), rA0
            #elif KB == 3
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 3
               movddup rB1, rB0
            #elif KB == 3
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 3
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd -16(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd -16(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 4
               movapd 0(pB), rB1
            #elif KB == 4
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 4
               movapd 0(pA), rA0
            #elif KB == 4
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 4
               movddup rB1, rB0
            #elif KB == 4
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 4
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 16(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 16(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 5
               movapd 32(pB), rB1
            #elif KB == 5
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 5
               movapd 32(pA), rA0
            #elif KB == 5
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 5
               movddup rB1, rB0
            #elif KB == 5
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 5
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 48(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 48(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 6
               movapd 64(pB), rB1
            #elif KB == 6
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 6
               movapd 64(pA), rA0
            #elif KB == 6
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 6
               movddup rB1, rB0
            #elif KB == 6
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 6
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 80(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 80(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 7
               movapd 96(pB), rB1
            #elif KB == 7
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 7
               movapd 96(pA), rA0
            #elif KB == 7
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 7
               movddup rB1, rB0
            #elif KB == 7
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 7
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 112(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 112(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 8
               movapd 128(pB), rB1
            #elif KB == 8
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 8
               movapd 128(pA), rA0
            #elif KB == 8
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 8
               movddup rB1, rB0
            #elif KB == 8
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 8
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 144(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 144(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 9
               movapd 160(pB), rB1
            #elif KB == 9
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 9
               movapd 160(pA), rA0
            #elif KB == 9
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 9
               movddup rB1, rB0
            #elif KB == 9
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 9
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 176(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 176(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 10
               movapd 192(pB), rB1
            #elif KB == 10
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 10
               movapd 192(pA), rA0
            #elif KB == 10
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 10
               movddup rB1, rB0
            #elif KB == 10
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 10
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 208(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 208(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 11
               movapd 224(pB), rB1
            #elif KB == 11
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 11
               movapd 224(pA), rA0
            #elif KB == 11
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 11
               movddup rB1, rB0
            #elif KB == 11
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 11
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 240(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 240(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 12
               movapd 256(pB), rB1
            #elif KB == 12
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 12
               movapd 256(pA), rA0
            #elif KB == 12
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 12
               movddup rB1, rB0
            #elif KB == 12
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 12
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 272(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 272(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 13
               movapd 288(pB), rB1
            #elif KB == 13
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 13
               movapd 288(pA), rA0
            #elif KB == 13
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 13
               movddup rB1, rB0
            #elif KB == 13
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 13
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 304(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 304(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 14
               movapd 320(pB), rB1
            #elif KB == 14
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 14
               movapd 320(pA), rA0
            #elif KB == 14
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 14
               movddup rB1, rB0
            #elif KB == 14
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 14
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 336(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 336(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 15
               movapd 352(pB), rB1
            #elif KB == 15
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 15
               movapd 352(pA), rA0
            #elif KB == 15
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 15
               movddup rB1, rB0
            #elif KB == 15
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 15
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 368(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 368(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 16
               movapd 384(pB), rB1
            #elif KB == 16
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 16
               movapd 384(pA), rA0
            #elif KB == 16
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 16
               movddup rB1, rB0
            #elif KB == 16
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 16
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 400(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 400(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 17
               movapd 416(pB), rB1
            #elif KB == 17
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 17
               movapd 416(pA), rA0
            #elif KB == 17
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 17
               movddup rB1, rB0
            #elif KB == 17
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 17
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 432(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 432(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 18
               movapd 448(pB), rB1
            #elif KB == 18
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 18
               movapd 448(pA), rA0
            #elif KB == 18
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 18
               movddup rB1, rB0
            #elif KB == 18
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 18
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 464(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 464(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 19
               movapd 480(pB), rB1
            #elif KB == 19
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 19
               movapd 480(pA), rA0
            #elif KB == 19
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 19
               movddup rB1, rB0
            #elif KB == 19
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 19
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 496(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 496(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 20
               movapd 512(pB), rB1
            #elif KB == 20
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 20
               movapd 512(pA), rA0
            #elif KB == 20
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 20
               movddup rB1, rB0
            #elif KB == 20
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 20
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 528(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 528(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 21
               movapd 544(pB), rB1
            #elif KB == 21
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 21
               movapd 544(pA), rA0
            #elif KB == 21
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 21
               movddup rB1, rB0
            #elif KB == 21
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 21
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 560(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 560(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 22
               movapd 576(pB), rB1
            #elif KB == 22
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 22
               movapd 576(pA), rA0
            #elif KB == 22
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 22
               movddup rB1, rB0
            #elif KB == 22
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 22
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 592(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 592(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 23
               movapd 608(pB), rB1
            #elif KB == 23
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 23
               movapd 608(pA), rA0
            #elif KB == 23
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 23
               movddup rB1, rB0
            #elif KB == 23
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 23
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 624(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 624(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 24
               movapd 640(pB), rB1
            #elif KB == 24
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 24
               movapd 640(pA), rA0
            #elif KB == 24
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 24
               movddup rB1, rB0
            #elif KB == 24
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 24
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 656(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 656(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 25
               movapd 672(pB), rB1
            #elif KB == 25
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 25
               movapd 672(pA), rA0
            #elif KB == 25
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 25
               movddup rB1, rB0
            #elif KB == 25
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 25
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 688(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 688(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 26
               movapd 704(pB), rB1
            #elif KB == 26
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 26
               movapd 704(pA), rA0
            #elif KB == 26
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 26
               movddup rB1, rB0
            #elif KB == 26
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 26
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 720(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 720(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 27
               movapd 736(pB), rB1
            #elif KB == 27
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 27
               movapd 736(pA), rA0
            #elif KB == 27
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 27
               movddup rB1, rB0
            #elif KB == 27
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 27
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 752(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 752(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 28
               movapd 768(pB), rB1
            #elif KB == 28
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 28
               movapd 768(pA), rA0
            #elif KB == 28
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 28
               movddup rB1, rB0
            #elif KB == 28
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 28
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 784(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 784(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 29
               movapd 800(pB), rB1
            #elif KB == 29
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 29
               movapd 800(pA), rA0
            #elif KB == 29
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 29
               movddup rB1, rB0
            #elif KB == 29
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 29
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 816(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 816(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 30
               movapd 832(pB), rB1
            #elif KB == 30
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 30
               movapd 832(pA), rA0
            #elif KB == 30
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 30
               movddup rB1, rB0
            #elif KB == 30
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 30
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 848(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 848(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 31
               movapd 864(pB), rB1
            #elif KB == 31
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 31
               movapd 864(pA), rA0
            #elif KB == 31
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 31
               movddup rB1, rB0
            #elif KB == 31
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 31
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 880(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 880(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 32
               movapd 896(pB), rB1
            #elif KB == 32
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 32
               movapd 896(pA), rA0
            #elif KB == 32
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 32
               movddup rB1, rB0
            #elif KB == 32
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 32
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 912(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 912(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 33
               movapd 928(pB), rB1
            #elif KB == 33
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 33
               movapd 928(pA), rA0
            #elif KB == 33
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 33
               movddup rB1, rB0
            #elif KB == 33
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 33
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 944(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 944(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 34
               movapd 960(pB), rB1
            #elif KB == 34
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 34
               movapd 960(pA), rA0
            #elif KB == 34
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 34
               movddup rB1, rB0
            #elif KB == 34
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 34
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 976(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 976(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 35
               movapd 992(pB), rB1
            #elif KB == 35
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 35
               movapd 992(pA), rA0
            #elif KB == 35
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 35
               movddup rB1, rB0
            #elif KB == 35
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 35
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1008(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1008(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 36
               movapd 1024(pB), rB1
            #elif KB == 36
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 36
               movapd 1024(pA), rA0
            #elif KB == 36
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 36
               movddup rB1, rB0
            #elif KB == 36
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 36
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1040(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1040(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 37
               movapd 1056(pB), rB1
            #elif KB == 37
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 37
               movapd 1056(pA), rA0
            #elif KB == 37
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 37
               movddup rB1, rB0
            #elif KB == 37
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 37
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1072(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1072(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 38
               movapd 1088(pB), rB1
            #elif KB == 38
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 38
               movapd 1088(pA), rA0
            #elif KB == 38
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 38
               movddup rB1, rB0
            #elif KB == 38
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 38
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1104(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1104(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 39
               movapd 1120(pB), rB1
            #elif KB == 39
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 39
               movapd 1120(pA), rA0
            #elif KB == 39
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 39
               movddup rB1, rB0
            #elif KB == 39
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 39
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1136(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1136(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 40
               movapd 1152(pB), rB1
            #elif KB == 40
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 40
               movapd 1152(pA), rA0
            #elif KB == 40
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 40
               movddup rB1, rB0
            #elif KB == 40
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 40
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1168(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1168(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 41
               movapd 1184(pB), rB1
            #elif KB == 41
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 41
               movapd 1184(pA), rA0
            #elif KB == 41
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 41
               movddup rB1, rB0
            #elif KB == 41
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 41
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1200(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1200(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 42
               movapd 1216(pB), rB1
            #elif KB == 42
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 42
               movapd 1216(pA), rA0
            #elif KB == 42
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 42
               movddup rB1, rB0
            #elif KB == 42
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 42
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1232(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1232(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 43
               movapd 1248(pB), rB1
            #elif KB == 43
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 43
               movapd 1248(pA), rA0
            #elif KB == 43
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 43
               movddup rB1, rB0
            #elif KB == 43
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 43
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1264(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1264(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 44
               movapd 1280(pB), rB1
            #elif KB == 44
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 44
               movapd 1280(pA), rA0
            #elif KB == 44
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 44
               movddup rB1, rB0
            #elif KB == 44
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 44
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1296(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1296(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 45
               movapd 1312(pB), rB1
            #elif KB == 45
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 45
               movapd 1312(pA), rA0
            #elif KB == 45
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 45
               movddup rB1, rB0
            #elif KB == 45
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 45
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1328(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1328(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 46
               movapd 1344(pB), rB1
            #elif KB == 46
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 46
               movapd 1344(pA), rA0
            #elif KB == 46
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 46
               movddup rB1, rB0
            #elif KB == 46
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 46
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1360(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1360(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 47
               movapd 1376(pB), rB1
            #elif KB == 47
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 47
               movapd 1376(pA), rA0
            #elif KB == 47
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 47
               movddup rB1, rB0
            #elif KB == 47
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 47
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1392(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1392(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 48
               movapd 1408(pB), rB1
            #elif KB == 48
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 48
               movapd 1408(pA), rA0
            #elif KB == 48
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 48
               movddup rB1, rB0
            #elif KB == 48
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 48
            nop
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1424(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1424(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 49
               movapd 1440(pB), rB1
            #elif KB == 49
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 49
               movapd 1440(pA), rA0
            #elif KB == 49
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 49
               movddup rB1, rB0
            #elif KB == 49
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 49
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1456(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1456(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 50
               movapd 1472(pB), rB1
            #elif KB == 50
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 50
               movapd 1472(pA), rA0
            #elif KB == 50
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 50
               movddup rB1, rB0
            #elif KB == 50
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 50
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1488(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1488(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 51
               movapd 1504(pB), rB1
            #elif KB == 51
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 51
               movapd 1504(pA), rA0
            #elif KB == 51
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 51
               movddup rB1, rB0
            #elif KB == 51
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 51
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1520(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1520(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 52
               movapd 1536(pB), rB1
            #elif KB == 52
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 52
               movapd 1536(pA), rA0
            #elif KB == 52
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 52
               movddup rB1, rB0
            #elif KB == 52
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 52
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1552(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1552(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 53
               movapd 1568(pB), rB1
            #elif KB == 53
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 53
               movapd 1568(pA), rA0
            #elif KB == 53
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 53
               movddup rB1, rB0
            #elif KB == 53
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 53
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1584(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1584(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 54
               movapd 1600(pB), rB1
            #elif KB == 54
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 54
               movapd 1600(pA), rA0
            #elif KB == 54
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 54
               movddup rB1, rB0
            #elif KB == 54
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 54
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1616(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1616(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 55
               movapd 1632(pB), rB1
            #elif KB == 55
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 55
               movapd 1632(pA), rA0
            #elif KB == 55
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 55
               movddup rB1, rB0
            #elif KB == 55
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 55
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1648(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1648(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 56
               movapd 1664(pB), rB1
            #elif KB == 56
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 56
               movapd 1664(pA), rA0
            #elif KB == 56
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 56
               movddup rB1, rB0
            #elif KB == 56
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 56
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1680(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1680(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 57
               movapd 1696(pB), rB1
            #elif KB == 57
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 57
               movapd 1696(pA), rA0
            #elif KB == 57
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 57
               movddup rB1, rB0
            #elif KB == 57
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 57
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1712(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1712(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 58
               movapd 1728(pB), rB1
            #elif KB == 58
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 58
               movapd 1728(pA), rA0
            #elif KB == 58
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 58
               movddup rB1, rB0
            #elif KB == 58
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 58
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1744(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1744(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 59
               movapd 1760(pB), rB1
            #elif KB == 59
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 59
               movapd 1760(pA), rA0
            #elif KB == 59
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 59
               movddup rB1, rB0
            #elif KB == 59
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 59
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1776(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1776(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 60
               movapd 1792(pB), rB1
            #elif KB == 60
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 60
               movapd 1792(pA), rA0
            #elif KB == 60
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 60
               movddup rB1, rB0
            #elif KB == 60
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 60
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1808(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1808(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 61
               movapd 1824(pB), rB1
            #elif KB == 61
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 61
               movapd 1824(pA), rA0
            #elif KB == 61
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 61
               movddup rB1, rB0
            #elif KB == 61
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 61
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1840(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1840(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 62
               movapd 1856(pB), rB1
            #elif KB == 62
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 62
               movapd 1856(pA), rA0
            #elif KB == 62
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 62
               movddup rB1, rB0
            #elif KB == 62
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 62
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1872(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1872(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 63
               movapd 1888(pB), rB1
            #elif KB == 63
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 63
               movapd 1888(pA), rA0
            #elif KB == 63
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 63
               movddup rB1, rB0
            #elif KB == 63
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 63
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1904(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1904(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 64
               movapd 1920(pB), rB1
            #elif KB == 64
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 64
               movapd 1920(pA), rA0
            #elif KB == 64
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 64
               movddup rB1, rB0
            #elif KB == 64
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 64
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1936(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1936(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 65
               movapd 1952(pB), rB1
            #elif KB == 65
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 65
               movapd 1952(pA), rA0
            #elif KB == 65
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 65
               movddup rB1, rB0
            #elif KB == 65
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 65
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 1968(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 1968(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 66
               movapd 1984(pB), rB1
            #elif KB == 66
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 66
               movapd 1984(pA), rA0
            #elif KB == 66
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 66
               movddup rB1, rB0
            #elif KB == 66
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 66
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2000(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2000(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 67
               movapd 2016(pB), rB1
            #elif KB == 67
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 67
               movapd 2016(pA), rA0
            #elif KB == 67
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 67
               movddup rB1, rB0
            #elif KB == 67
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 67
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2032(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2032(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 68
               movapd 2048(pB), rB1
            #elif KB == 68
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 68
               movapd 2048(pA), rA0
            #elif KB == 68
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 68
               movddup rB1, rB0
            #elif KB == 68
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 68
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2064(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2064(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 69
               movapd 2080(pB), rB1
            #elif KB == 69
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 69
               movapd 2080(pA), rA0
            #elif KB == 69
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 69
               movddup rB1, rB0
            #elif KB == 69
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 69
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2096(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2096(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 70
               movapd 2112(pB), rB1
            #elif KB == 70
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 70
               movapd 2112(pA), rA0
            #elif KB == 70
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 70
               movddup rB1, rB0
            #elif KB == 70
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 70
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2128(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2128(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 71
               movapd 2144(pB), rB1
            #elif KB == 71
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 71
               movapd 2144(pA), rA0
            #elif KB == 71
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 71
               movddup rB1, rB0
            #elif KB == 71
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 71
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2160(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2160(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 72
               movapd 2176(pB), rB1
            #elif KB == 72
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 72
               movapd 2176(pA), rA0
            #elif KB == 72
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 72
               movddup rB1, rB0
            #elif KB == 72
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 72
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2192(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2192(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 73
               movapd 2208(pB), rB1
            #elif KB == 73
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 73
               movapd 2208(pA), rA0
            #elif KB == 73
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 73
               movddup rB1, rB0
            #elif KB == 73
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 73
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2224(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2224(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 74
               movapd 2240(pB), rB1
            #elif KB == 74
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 74
               movapd 2240(pA), rA0
            #elif KB == 74
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 74
               movddup rB1, rB0
            #elif KB == 74
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 74
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2256(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2256(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 75
               movapd 2272(pB), rB1
            #elif KB == 75
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 75
               movapd 2272(pA), rA0
            #elif KB == 75
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 75
               movddup rB1, rB0
            #elif KB == 75
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 75
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2288(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2288(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 76
               movapd 2304(pB), rB1
            #elif KB == 76
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 76
               movapd 2304(pA), rA0
            #elif KB == 76
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 76
               movddup rB1, rB0
            #elif KB == 76
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 76
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2320(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2320(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 77
               movapd 2336(pB), rB1
            #elif KB == 77
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 77
               movapd 2336(pA), rA0
            #elif KB == 77
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 77
               movddup rB1, rB0
            #elif KB == 77
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 77
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2352(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2352(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 78
               movapd 2368(pB), rB1
            #elif KB == 78
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 78
               movapd 2368(pA), rA0
            #elif KB == 78
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 78
               movddup rB1, rB0
            #elif KB == 78
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 78
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2384(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2384(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 79
               movapd 2400(pB), rB1
            #elif KB == 79
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 79
               movapd 2400(pA), rA0
            #elif KB == 79
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 79
               movddup rB1, rB0
            #elif KB == 79
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 79
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2416(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2416(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 80
               movapd 2432(pB), rB1
            #elif KB == 80
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 80
               movapd 2432(pA), rA0
            #elif KB == 80
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 80
               movddup rB1, rB0
            #elif KB == 80
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 80
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2448(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2448(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 81
               movapd 2464(pB), rB1
            #elif KB == 81
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 81
               movapd 2464(pA), rA0
            #elif KB == 81
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 81
               movddup rB1, rB0
            #elif KB == 81
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 81
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2480(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2480(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 82
               movapd 2496(pB), rB1
            #elif KB == 82
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 82
               movapd 2496(pA), rA0
            #elif KB == 82
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 82
               movddup rB1, rB0
            #elif KB == 82
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 82
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2512(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2512(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 83
               movapd 2528(pB), rB1
            #elif KB == 83
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 83
               movapd 2528(pA), rA0
            #elif KB == 83
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 83
               movddup rB1, rB0
            #elif KB == 83
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 83
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2544(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2544(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 84
               movapd 2560(pB), rB1
            #elif KB == 84
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 84
               movapd 2560(pA), rA0
            #elif KB == 84
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 84
               movddup rB1, rB0
            #elif KB == 84
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 84
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2576(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2576(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 85
               movapd 2592(pB), rB1
            #elif KB == 85
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 85
               movapd 2592(pA), rA0
            #elif KB == 85
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 85
               movddup rB1, rB0
            #elif KB == 85
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 85
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2608(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2608(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 86
               movapd 2624(pB), rB1
            #elif KB == 86
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 86
               movapd 2624(pA), rA0
            #elif KB == 86
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 86
               movddup rB1, rB0
            #elif KB == 86
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 86
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2640(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2640(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 87
               movapd 2656(pB), rB1
            #elif KB == 87
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 87
               movapd 2656(pA), rA0
            #elif KB == 87
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 87
               movddup rB1, rB0
            #elif KB == 87
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 87
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2672(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2672(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 88
               movapd 2688(pB), rB1
            #elif KB == 88
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 88
               movapd 2688(pA), rA0
            #elif KB == 88
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 88
               movddup rB1, rB0
            #elif KB == 88
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 88
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2704(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2704(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 89
               movapd 2720(pB), rB1
            #elif KB == 89
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 89
               movapd 2720(pA), rA0
            #elif KB == 89
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 89
               movddup rB1, rB0
            #elif KB == 89
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 89
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2736(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2736(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 90
               movapd 2752(pB), rB1
            #elif KB == 90
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 90
               movapd 2752(pA), rA0
            #elif KB == 90
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 90
               movddup rB1, rB0
            #elif KB == 90
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 90
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2768(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2768(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 91
               movapd 2784(pB), rB1
            #elif KB == 91
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 91
               movapd 2784(pA), rA0
            #elif KB == 91
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 91
               movddup rB1, rB0
            #elif KB == 91
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 91
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2800(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2800(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 92
               movapd 2816(pB), rB1
            #elif KB == 92
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 92
               movapd 2816(pA), rA0
            #elif KB == 92
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 92
               movddup rB1, rB0
            #elif KB == 92
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 92
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2832(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2832(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 93
               movapd 2848(pB), rB1
            #elif KB == 93
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 93
               movapd 2848(pA), rA0
            #elif KB == 93
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 93
               movddup rB1, rB0
            #elif KB == 93
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 93
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2864(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2864(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 94
               movapd 2880(pB), rB1
            #elif KB == 94
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 94
               movapd 2880(pA), rA0
            #elif KB == 94
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 94
               movddup rB1, rB0
            #elif KB == 94
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 94
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2896(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2896(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 95
               movapd 2912(pB), rB1
            #elif KB == 95
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 95
               movapd 2912(pA), rA0
            #elif KB == 95
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 95
               movddup rB1, rB0
            #elif KB == 95
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 95
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2928(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2928(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 96
               movapd 2944(pB), rB1
            #elif KB == 96
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 96
               movapd 2944(pA), rA0
            #elif KB == 96
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 96
               movddup rB1, rB0
            #elif KB == 96
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 96
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2960(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2960(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 97
               movapd 2976(pB), rB1
            #elif KB == 97
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 97
               movapd 2976(pA), rA0
            #elif KB == 97
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 97
               movddup rB1, rB0
            #elif KB == 97
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 97
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 2992(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 2992(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 98
               movapd 3008(pB), rB1
            #elif KB == 98
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 98
               movapd 3008(pA), rA0
            #elif KB == 98
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 98
               movddup rB1, rB0
            #elif KB == 98
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 98
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3024(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3024(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 99
               movapd 3040(pB), rB1
            #elif KB == 99
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 99
               movapd 3040(pA), rA0
            #elif KB == 99
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 99
               movddup rB1, rB0
            #elif KB == 99
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 99
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3056(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3056(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 100
               movapd 3072(pB), rB1
            #elif KB == 100
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 100
               movapd 3072(pA), rA0
            #elif KB == 100
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 100
               movddup rB1, rB0
            #elif KB == 100
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 100
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3088(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3088(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 101
               movapd 3104(pB), rB1
            #elif KB == 101
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 101
               movapd 3104(pA), rA0
            #elif KB == 101
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 101
               movddup rB1, rB0
            #elif KB == 101
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 101
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3120(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3120(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 102
               movapd 3136(pB), rB1
            #elif KB == 102
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 102
               movapd 3136(pA), rA0
            #elif KB == 102
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 102
               movddup rB1, rB0
            #elif KB == 102
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 102
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3152(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3152(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 103
               movapd 3168(pB), rB1
            #elif KB == 103
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 103
               movapd 3168(pA), rA0
            #elif KB == 103
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 103
               movddup rB1, rB0
            #elif KB == 103
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 103
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3184(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3184(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 104
               movapd 3200(pB), rB1
            #elif KB == 104
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 104
               movapd 3200(pA), rA0
            #elif KB == 104
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 104
               movddup rB1, rB0
            #elif KB == 104
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 104
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3216(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3216(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 105
               movapd 3232(pB), rB1
            #elif KB == 105
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 105
               movapd 3232(pA), rA0
            #elif KB == 105
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 105
               movddup rB1, rB0
            #elif KB == 105
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 105
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3248(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3248(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 106
               movapd 3264(pB), rB1
            #elif KB == 106
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 106
               movapd 3264(pA), rA0
            #elif KB == 106
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 106
               movddup rB1, rB0
            #elif KB == 106
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 106
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3280(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3280(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 107
               movapd 3296(pB), rB1
            #elif KB == 107
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 107
               movapd 3296(pA), rA0
            #elif KB == 107
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 107
               movddup rB1, rB0
            #elif KB == 107
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 107
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3312(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3312(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 108
               movapd 3328(pB), rB1
            #elif KB == 108
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 108
               movapd 3328(pA), rA0
            #elif KB == 108
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 108
               movddup rB1, rB0
            #elif KB == 108
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 108
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3344(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3344(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 109
               movapd 3360(pB), rB1
            #elif KB == 109
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 109
               movapd 3360(pA), rA0
            #elif KB == 109
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 109
               movddup rB1, rB0
            #elif KB == 109
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 109
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3376(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3376(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 110
               movapd 3392(pB), rB1
            #elif KB == 110
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 110
               movapd 3392(pA), rA0
            #elif KB == 110
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 110
               movddup rB1, rB0
            #elif KB == 110
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 110
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3408(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3408(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 111
               movapd 3424(pB), rB1
            #elif KB == 111
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 111
               movapd 3424(pA), rA0
            #elif KB == 111
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 111
               movddup rB1, rB0
            #elif KB == 111
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 111
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3440(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3440(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 112
               movapd 3456(pB), rB1
            #elif KB == 112
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 112
               movapd 3456(pA), rA0
            #elif KB == 112
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 112
               movddup rB1, rB0
            #elif KB == 112
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 112
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3472(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3472(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 113
               movapd 3488(pB), rB1
            #elif KB == 113
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 113
               movapd 3488(pA), rA0
            #elif KB == 113
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 113
               movddup rB1, rB0
            #elif KB == 113
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 113
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3504(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3504(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 114
               movapd 3520(pB), rB1
            #elif KB == 114
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 114
               movapd 3520(pA), rA0
            #elif KB == 114
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 114
               movddup rB1, rB0
            #elif KB == 114
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 114
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3536(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3536(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 115
               movapd 3552(pB), rB1
            #elif KB == 115
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 115
               movapd 3552(pA), rA0
            #elif KB == 115
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 115
               movddup rB1, rB0
            #elif KB == 115
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 115
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3568(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3568(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 116
               movapd 3584(pB), rB1
            #elif KB == 116
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 116
               movapd 3584(pA), rA0
            #elif KB == 116
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 116
               movddup rB1, rB0
            #elif KB == 116
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 116
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3600(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3600(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 117
               movapd 3616(pB), rB1
            #elif KB == 117
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 117
               movapd 3616(pA), rA0
            #elif KB == 117
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 117
               movddup rB1, rB0
            #elif KB == 117
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 117
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3632(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3632(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 118
               movapd 3648(pB), rB1
            #elif KB == 118
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 118
               movapd 3648(pA), rA0
            #elif KB == 118
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 118
               movddup rB1, rB0
            #elif KB == 118
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 118
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3664(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3664(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 119
               movapd 3680(pB), rB1
            #elif KB == 119
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 119
               movapd 3680(pA), rA0
            #elif KB == 119
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 119
               movddup rB1, rB0
            #elif KB == 119
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 119
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3696(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3696(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 120
               movapd 3712(pB), rB1
            #elif KB == 120
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 120
               movapd 3712(pA), rA0
            #elif KB == 120
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 120
               movddup rB1, rB0
            #elif KB == 120
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 120
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3728(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3728(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 121
               movapd 3744(pB), rB1
            #elif KB == 121
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 121
               movapd 3744(pA), rA0
            #elif KB == 121
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 121
               movddup rB1, rB0
            #elif KB == 121
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 121
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3760(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3760(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 122
               movapd 3776(pB), rB1
            #elif KB == 122
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 122
               movapd 3776(pA), rA0
            #elif KB == 122
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 122
               movddup rB1, rB0
            #elif KB == 122
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 122
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3792(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3792(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 123
               movapd 3808(pB), rB1
            #elif KB == 123
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 123
               movapd 3808(pA), rA0
            #elif KB == 123
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 123
               movddup rB1, rB0
            #elif KB == 123
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 123
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3824(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3824(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 124
               movapd 3840(pB), rB1
            #elif KB == 124
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 124
               movapd 3840(pA), rA0
            #elif KB == 124
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 124
               movddup rB1, rB0
            #elif KB == 124
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 124
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3856(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3856(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 125
               movapd 3872(pB), rB1
            #elif KB == 125
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 125
               movapd 3872(pA), rA0
            #elif KB == 125
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 125
               movddup rB1, rB0
            #elif KB == 125
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 125
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3888(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3888(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 126
               movapd 3904(pB), rB1
            #elif KB == 126
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 126
               movapd 3904(pA), rA0
            #elif KB == 126
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 126
               movddup rB1, rB0
            #elif KB == 126
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 126
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3920(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3920(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 127
               movapd 3936(pB), rB1
            #elif KB == 127
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 127
               movapd 3936(pA), rA0
            #elif KB == 127
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 127
               movddup rB1, rB0
            #elif KB == 127
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 127
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3952(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3952(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 128
               movapd 3968(pB), rB1
            #elif KB == 128
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 128
               movapd 3968(pA), rA0
            #elif KB == 128
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 128
               movddup rB1, rB0
            #elif KB == 128
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 128
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 3984(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 3984(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 129
               movapd 4000(pB), rB1
            #elif KB == 129
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 129
               movapd 4000(pA), rA0
            #elif KB == 129
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 129
               movddup rB1, rB0
            #elif KB == 129
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 129
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4016(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4016(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 130
               movapd 4032(pB), rB1
            #elif KB == 130
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 130
               movapd 4032(pA), rA0
            #elif KB == 130
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 130
               movddup rB1, rB0
            #elif KB == 130
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 130
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4048(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4048(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 131
               movapd 4064(pB), rB1
            #elif KB == 131
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 131
               movapd 4064(pA), rA0
            #elif KB == 131
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 131
               movddup rB1, rB0
            #elif KB == 131
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 131
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4080(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4080(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 132
               movapd 4096(pB), rB1
            #elif KB == 132
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 132
               movapd 4096(pA), rA0
            #elif KB == 132
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 132
               movddup rB1, rB0
            #elif KB == 132
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 132
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4112(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4112(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 133
               movapd 4128(pB), rB1
            #elif KB == 133
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 133
               movapd 4128(pA), rA0
            #elif KB == 133
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 133
               movddup rB1, rB0
            #elif KB == 133
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 133
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4144(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4144(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 134
               movapd 4160(pB), rB1
            #elif KB == 134
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 134
               movapd 4160(pA), rA0
            #elif KB == 134
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 134
               movddup rB1, rB0
            #elif KB == 134
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 134
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4176(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4176(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 135
               movapd 4192(pB), rB1
            #elif KB == 135
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 135
               movapd 4192(pA), rA0
            #elif KB == 135
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 135
               movddup rB1, rB0
            #elif KB == 135
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 135
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4208(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4208(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 136
               movapd 4224(pB), rB1
            #elif KB == 136
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 136
               movapd 4224(pA), rA0
            #elif KB == 136
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 136
               movddup rB1, rB0
            #elif KB == 136
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 136
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4240(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4240(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 137
               movapd 4256(pB), rB1
            #elif KB == 137
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 137
               movapd 4256(pA), rA0
            #elif KB == 137
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 137
               movddup rB1, rB0
            #elif KB == 137
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 137
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4272(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4272(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 138
               movapd 4288(pB), rB1
            #elif KB == 138
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 138
               movapd 4288(pA), rA0
            #elif KB == 138
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 138
               movddup rB1, rB0
            #elif KB == 138
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 138
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4304(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4304(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 139
               movapd 4320(pB), rB1
            #elif KB == 139
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 139
               movapd 4320(pA), rA0
            #elif KB == 139
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 139
               movddup rB1, rB0
            #elif KB == 139
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 139
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4336(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4336(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 140
               movapd 4352(pB), rB1
            #elif KB == 140
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 140
               movapd 4352(pA), rA0
            #elif KB == 140
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 140
               movddup rB1, rB0
            #elif KB == 140
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 140
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4368(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4368(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 141
               movapd 4384(pB), rB1
            #elif KB == 141
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 141
               movapd 4384(pA), rA0
            #elif KB == 141
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 141
               movddup rB1, rB0
            #elif KB == 141
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 141
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4400(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4400(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 142
               movapd 4416(pB), rB1
            #elif KB == 142
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 142
               movapd 4416(pA), rA0
            #elif KB == 142
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 142
               movddup rB1, rB0
            #elif KB == 142
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 142
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4432(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4432(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 143
               movapd 4448(pB), rB1
            #elif KB == 143
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 143
               movapd 4448(pA), rA0
            #elif KB == 143
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 143
               movddup rB1, rB0
            #elif KB == 143
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 143
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4464(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4464(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 144
               movapd 4480(pB), rB1
            #elif KB == 144
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 144
               movapd 4480(pA), rA0
            #elif KB == 144
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 144
               movddup rB1, rB0
            #elif KB == 144
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 144
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4496(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4496(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 145
               movapd 4512(pB), rB1
            #elif KB == 145
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 145
               movapd 4512(pA), rA0
            #elif KB == 145
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 145
               movddup rB1, rB0
            #elif KB == 145
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 145
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4528(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4528(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 146
               movapd 4544(pB), rB1
            #elif KB == 146
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 146
               movapd 4544(pA), rA0
            #elif KB == 146
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 146
               movddup rB1, rB0
            #elif KB == 146
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 146
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4560(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4560(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 147
               movapd 4576(pB), rB1
            #elif KB == 147
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 147
               movapd 4576(pA), rA0
            #elif KB == 147
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 147
               movddup rB1, rB0
            #elif KB == 147
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 147
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4592(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4592(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 148
               movapd 4608(pB), rB1
            #elif KB == 148
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 148
               movapd 4608(pA), rA0
            #elif KB == 148
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 148
               movddup rB1, rB0
            #elif KB == 148
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 148
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4624(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4624(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 149
               movapd 4640(pB), rB1
            #elif KB == 149
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 149
               movapd 4640(pA), rA0
            #elif KB == 149
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 149
               movddup rB1, rB0
            #elif KB == 149
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 149
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4656(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4656(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 150
               movapd 4672(pB), rB1
            #elif KB == 150
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 150
               movapd 4672(pA), rA0
            #elif KB == 150
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 150
               movddup rB1, rB0
            #elif KB == 150
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 150
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4688(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4688(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 151
               movapd 4704(pB), rB1
            #elif KB == 151
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 151
               movapd 4704(pA), rA0
            #elif KB == 151
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 151
               movddup rB1, rB0
            #elif KB == 151
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 151
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4720(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4720(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 152
               movapd 4736(pB), rB1
            #elif KB == 152
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 152
               movapd 4736(pA), rA0
            #elif KB == 152
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 152
               movddup rB1, rB0
            #elif KB == 152
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 152
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4752(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4752(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 153
               movapd 4768(pB), rB1
            #elif KB == 153
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 153
               movapd 4768(pA), rA0
            #elif KB == 153
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 153
               movddup rB1, rB0
            #elif KB == 153
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 153
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4784(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4784(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 154
               movapd 4800(pB), rB1
            #elif KB == 154
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 154
               movapd 4800(pA), rA0
            #elif KB == 154
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 154
               movddup rB1, rB0
            #elif KB == 154
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 154
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4816(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4816(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 155
               movapd 4832(pB), rB1
            #elif KB == 155
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 155
               movapd 4832(pA), rA0
            #elif KB == 155
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 155
               movddup rB1, rB0
            #elif KB == 155
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 155
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4848(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4848(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 156
               movapd 4864(pB), rB1
            #elif KB == 156
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 156
               movapd 4864(pA), rA0
            #elif KB == 156
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 156
               movddup rB1, rB0
            #elif KB == 156
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 156
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4880(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4880(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 157
               movapd 4896(pB), rB1
            #elif KB == 157
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 157
               movapd 4896(pA), rA0
            #elif KB == 157
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 157
               movddup rB1, rB0
            #elif KB == 157
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 157
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4912(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4912(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 158
               movapd 4928(pB), rB1
            #elif KB == 158
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 158
               movapd 4928(pA), rA0
            #elif KB == 158
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 158
               movddup rB1, rB0
            #elif KB == 158
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 158
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4944(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4944(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 159
               movapd 4960(pB), rB1
            #elif KB == 159
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 159
               movapd 4960(pA), rA0
            #elif KB == 159
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 159
               movddup rB1, rB0
            #elif KB == 159
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 159
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 4976(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 4976(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 160
               movapd 4992(pB), rB1
            #elif KB == 160
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 160
               movapd 4992(pA), rA0
            #elif KB == 160
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 160
               movddup rB1, rB0
            #elif KB == 160
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 160
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5008(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5008(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 161
               movapd 5024(pB), rB1
            #elif KB == 161
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 161
               movapd 5024(pA), rA0
            #elif KB == 161
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 161
               movddup rB1, rB0
            #elif KB == 161
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 161
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5040(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5040(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 162
               movapd 5056(pB), rB1
            #elif KB == 162
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 162
               movapd 5056(pA), rA0
            #elif KB == 162
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 162
               movddup rB1, rB0
            #elif KB == 162
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 162
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5072(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5072(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 163
               movapd 5088(pB), rB1
            #elif KB == 163
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 163
               movapd 5088(pA), rA0
            #elif KB == 163
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 163
               movddup rB1, rB0
            #elif KB == 163
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 163
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5104(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5104(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 164
               movapd 5120(pB), rB1
            #elif KB == 164
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 164
               movapd 5120(pA), rA0
            #elif KB == 164
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 164
               movddup rB1, rB0
            #elif KB == 164
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 164
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5136(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5136(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 165
               movapd 5152(pB), rB1
            #elif KB == 165
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 165
               movapd 5152(pA), rA0
            #elif KB == 165
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 165
               movddup rB1, rB0
            #elif KB == 165
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 165
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5168(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5168(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 166
               movapd 5184(pB), rB1
            #elif KB == 166
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 166
               movapd 5184(pA), rA0
            #elif KB == 166
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 166
               movddup rB1, rB0
            #elif KB == 166
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 166
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5200(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5200(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 167
               movapd 5216(pB), rB1
            #elif KB == 167
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 167
               movapd 5216(pA), rA0
            #elif KB == 167
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 167
               movddup rB1, rB0
            #elif KB == 167
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 167
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5232(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5232(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 168
               movapd 5248(pB), rB1
            #elif KB == 168
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 168
               movapd 5248(pA), rA0
            #elif KB == 168
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 168
               movddup rB1, rB0
            #elif KB == 168
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 168
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5264(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5264(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 169
               movapd 5280(pB), rB1
            #elif KB == 169
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 169
               movapd 5280(pA), rA0
            #elif KB == 169
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 169
               movddup rB1, rB0
            #elif KB == 169
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 169
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5296(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5296(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 170
               movapd 5312(pB), rB1
            #elif KB == 170
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 170
               movapd 5312(pA), rA0
            #elif KB == 170
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 170
               movddup rB1, rB0
            #elif KB == 170
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 170
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5328(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5328(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 171
               movapd 5344(pB), rB1
            #elif KB == 171
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 171
               movapd 5344(pA), rA0
            #elif KB == 171
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 171
               movddup rB1, rB0
            #elif KB == 171
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 171
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5360(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5360(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 172
               movapd 5376(pB), rB1
            #elif KB == 172
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 172
               movapd 5376(pA), rA0
            #elif KB == 172
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 172
               movddup rB1, rB0
            #elif KB == 172
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 172
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5392(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5392(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 173
               movapd 5408(pB), rB1
            #elif KB == 173
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 173
               movapd 5408(pA), rA0
            #elif KB == 173
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 173
               movddup rB1, rB0
            #elif KB == 173
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 173
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5424(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5424(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 174
               movapd 5440(pB), rB1
            #elif KB == 174
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 174
               movapd 5440(pA), rA0
            #elif KB == 174
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 174
               movddup rB1, rB0
            #elif KB == 174
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 174
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5456(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5456(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 175
               movapd 5472(pB), rB1
            #elif KB == 175
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 175
               movapd 5472(pA), rA0
            #elif KB == 175
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 175
               movddup rB1, rB0
            #elif KB == 175
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 175
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5488(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5488(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 176
               movapd 5504(pB), rB1
            #elif KB == 176
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 176
               movapd 5504(pA), rA0
            #elif KB == 176
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 176
               movddup rB1, rB0
            #elif KB == 176
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 176
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5520(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5520(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 177
               movapd 5536(pB), rB1
            #elif KB == 177
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 177
               movapd 5536(pA), rA0
            #elif KB == 177
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 177
               movddup rB1, rB0
            #elif KB == 177
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 177
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5552(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5552(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 178
               movapd 5568(pB), rB1
            #elif KB == 178
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 178
               movapd 5568(pA), rA0
            #elif KB == 178
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 178
               movddup rB1, rB0
            #elif KB == 178
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 178
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5584(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5584(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 179
               movapd 5600(pB), rB1
            #elif KB == 179
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 179
               movapd 5600(pA), rA0
            #elif KB == 179
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 179
               movddup rB1, rB0
            #elif KB == 179
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 179
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5616(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5616(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 180
               movapd 5632(pB), rB1
            #elif KB == 180
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 180
               movapd 5632(pA), rA0
            #elif KB == 180
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 180
               movddup rB1, rB0
            #elif KB == 180
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 180
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5648(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5648(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 181
               movapd 5664(pB), rB1
            #elif KB == 181
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 181
               movapd 5664(pA), rA0
            #elif KB == 181
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 181
               movddup rB1, rB0
            #elif KB == 181
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 181
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5680(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5680(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 182
               movapd 5696(pB), rB1
            #elif KB == 182
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 182
               movapd 5696(pA), rA0
            #elif KB == 182
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 182
               movddup rB1, rB0
            #elif KB == 182
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 182
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5712(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5712(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 183
               movapd 5728(pB), rB1
            #elif KB == 183
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 183
               movapd 5728(pA), rA0
            #elif KB == 183
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 183
               movddup rB1, rB0
            #elif KB == 183
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 183
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5744(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5744(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 184
               movapd 5760(pB), rB1
            #elif KB == 184
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 184
               movapd 5760(pA), rA0
            #elif KB == 184
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 184
               movddup rB1, rB0
            #elif KB == 184
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 184
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5776(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5776(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 185
               movapd 5792(pB), rB1
            #elif KB == 185
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 185
               movapd 5792(pA), rA0
            #elif KB == 185
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 185
               movddup rB1, rB0
            #elif KB == 185
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 185
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5808(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5808(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 186
               movapd 5824(pB), rB1
            #elif KB == 186
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 186
               movapd 5824(pA), rA0
            #elif KB == 186
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 186
               movddup rB1, rB0
            #elif KB == 186
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 186
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5840(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5840(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 187
               movapd 5856(pB), rB1
            #elif KB == 187
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 187
               movapd 5856(pA), rA0
            #elif KB == 187
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 187
               movddup rB1, rB0
            #elif KB == 187
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 187
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5872(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5872(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 188
               movapd 5888(pB), rB1
            #elif KB == 188
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 188
               movapd 5888(pA), rA0
            #elif KB == 188
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 188
               movddup rB1, rB0
            #elif KB == 188
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 188
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5904(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5904(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 189
               movapd 5920(pB), rB1
            #elif KB == 189
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 189
               movapd 5920(pA), rA0
            #elif KB == 189
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 189
               movddup rB1, rB0
            #elif KB == 189
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 189
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5936(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5936(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 190
               movapd 5952(pB), rB1
            #elif KB == 190
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 190
               movapd 5952(pA), rA0
            #elif KB == 190
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 190
               movddup rB1, rB0
            #elif KB == 190
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 190
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 5968(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 5968(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 191
               movapd 5984(pB), rB1
            #elif KB == 191
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 191
               movapd 5984(pA), rA0
            #elif KB == 191
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 191
               movddup rB1, rB0
            #elif KB == 191
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 191
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6000(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6000(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 192
               movapd 6016(pB), rB1
            #elif KB == 192
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 192
               movapd 6016(pA), rA0
            #elif KB == 192
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 192
               movddup rB1, rB0
            #elif KB == 192
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 192
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6032(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6032(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 193
               movapd 6048(pB), rB1
            #elif KB == 193
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 193
               movapd 6048(pA), rA0
            #elif KB == 193
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 193
               movddup rB1, rB0
            #elif KB == 193
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 193
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6064(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6064(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 194
               movapd 6080(pB), rB1
            #elif KB == 194
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 194
               movapd 6080(pA), rA0
            #elif KB == 194
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 194
               movddup rB1, rB0
            #elif KB == 194
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 194
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6096(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6096(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 195
               movapd 6112(pB), rB1
            #elif KB == 195
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 195
               movapd 6112(pA), rA0
            #elif KB == 195
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 195
               movddup rB1, rB0
            #elif KB == 195
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 195
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6128(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6128(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 196
               movapd 6144(pB), rB1
            #elif KB == 196
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 196
               movapd 6144(pA), rA0
            #elif KB == 196
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 196
               movddup rB1, rB0
            #elif KB == 196
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 196
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6160(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6160(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 197
               movapd 6176(pB), rB1
            #elif KB == 197
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 197
               movapd 6176(pA), rA0
            #elif KB == 197
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 197
               movddup rB1, rB0
            #elif KB == 197
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 197
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6192(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6192(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 198
               movapd 6208(pB), rB1
            #elif KB == 198
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 198
               movapd 6208(pA), rA0
            #elif KB == 198
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 198
               movddup rB1, rB0
            #elif KB == 198
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 198
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6224(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6224(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 199
               movapd 6240(pB), rB1
            #elif KB == 199
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 199
               movapd 6240(pA), rA0
            #elif KB == 199
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 199
               movddup rB1, rB0
            #elif KB == 199
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 199
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6256(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6256(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 200
               movapd 6272(pB), rB1
            #elif KB == 200
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 200
               movapd 6272(pA), rA0
            #elif KB == 200
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 200
               movddup rB1, rB0
            #elif KB == 200
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 200
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6288(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6288(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 201
               movapd 6304(pB), rB1
            #elif KB == 201
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 201
               movapd 6304(pA), rA0
            #elif KB == 201
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 201
               movddup rB1, rB0
            #elif KB == 201
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 201
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6320(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6320(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 202
               movapd 6336(pB), rB1
            #elif KB == 202
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 202
               movapd 6336(pA), rA0
            #elif KB == 202
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 202
               movddup rB1, rB0
            #elif KB == 202
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 202
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6352(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6352(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 203
               movapd 6368(pB), rB1
            #elif KB == 203
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 203
               movapd 6368(pA), rA0
            #elif KB == 203
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 203
               movddup rB1, rB0
            #elif KB == 203
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 203
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6384(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6384(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 204
               movapd 6400(pB), rB1
            #elif KB == 204
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 204
               movapd 6400(pA), rA0
            #elif KB == 204
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 204
               movddup rB1, rB0
            #elif KB == 204
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 204
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6416(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6416(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 205
               movapd 6432(pB), rB1
            #elif KB == 205
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 205
               movapd 6432(pA), rA0
            #elif KB == 205
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 205
               movddup rB1, rB0
            #elif KB == 205
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 205
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6448(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6448(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 206
               movapd 6464(pB), rB1
            #elif KB == 206
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 206
               movapd 6464(pA), rA0
            #elif KB == 206
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 206
               movddup rB1, rB0
            #elif KB == 206
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 206
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6480(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6480(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 207
               movapd 6496(pB), rB1
            #elif KB == 207
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 207
               movapd 6496(pA), rA0
            #elif KB == 207
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 207
               movddup rB1, rB0
            #elif KB == 207
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 207
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6512(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6512(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 208
               movapd 6528(pB), rB1
            #elif KB == 208
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 208
               movapd 6528(pA), rA0
            #elif KB == 208
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 208
               movddup rB1, rB0
            #elif KB == 208
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 208
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6544(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6544(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 209
               movapd 6560(pB), rB1
            #elif KB == 209
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 209
               movapd 6560(pA), rA0
            #elif KB == 209
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 209
               movddup rB1, rB0
            #elif KB == 209
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 209
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6576(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6576(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 210
               movapd 6592(pB), rB1
            #elif KB == 210
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 210
               movapd 6592(pA), rA0
            #elif KB == 210
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 210
               movddup rB1, rB0
            #elif KB == 210
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 210
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6608(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6608(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 211
               movapd 6624(pB), rB1
            #elif KB == 211
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 211
               movapd 6624(pA), rA0
            #elif KB == 211
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 211
               movddup rB1, rB0
            #elif KB == 211
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 211
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6640(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6640(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 212
               movapd 6656(pB), rB1
            #elif KB == 212
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 212
               movapd 6656(pA), rA0
            #elif KB == 212
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 212
               movddup rB1, rB0
            #elif KB == 212
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 212
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6672(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6672(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 213
               movapd 6688(pB), rB1
            #elif KB == 213
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 213
               movapd 6688(pA), rA0
            #elif KB == 213
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 213
               movddup rB1, rB0
            #elif KB == 213
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 213
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6704(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6704(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 214
               movapd 6720(pB), rB1
            #elif KB == 214
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 214
               movapd 6720(pA), rA0
            #elif KB == 214
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 214
               movddup rB1, rB0
            #elif KB == 214
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 214
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6736(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6736(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 215
               movapd 6752(pB), rB1
            #elif KB == 215
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 215
               movapd 6752(pA), rA0
            #elif KB == 215
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 215
               movddup rB1, rB0
            #elif KB == 215
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 215
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6768(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6768(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 216
               movapd 6784(pB), rB1
            #elif KB == 216
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 216
               movapd 6784(pA), rA0
            #elif KB == 216
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 216
               movddup rB1, rB0
            #elif KB == 216
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 216
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6800(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6800(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 217
               movapd 6816(pB), rB1
            #elif KB == 217
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 217
               movapd 6816(pA), rA0
            #elif KB == 217
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 217
               movddup rB1, rB0
            #elif KB == 217
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 217
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6832(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6832(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 218
               movapd 6848(pB), rB1
            #elif KB == 218
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 218
               movapd 6848(pA), rA0
            #elif KB == 218
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 218
               movddup rB1, rB0
            #elif KB == 218
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 218
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6864(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6864(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 219
               movapd 6880(pB), rB1
            #elif KB == 219
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 219
               movapd 6880(pA), rA0
            #elif KB == 219
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 219
               movddup rB1, rB0
            #elif KB == 219
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 219
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6896(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6896(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 220
               movapd 6912(pB), rB1
            #elif KB == 220
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 220
               movapd 6912(pA), rA0
            #elif KB == 220
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 220
               movddup rB1, rB0
            #elif KB == 220
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 220
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6928(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6928(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 221
               movapd 6944(pB), rB1
            #elif KB == 221
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 221
               movapd 6944(pA), rA0
            #elif KB == 221
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 221
               movddup rB1, rB0
            #elif KB == 221
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 221
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6960(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6960(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 222
               movapd 6976(pB), rB1
            #elif KB == 222
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 222
               movapd 6976(pA), rA0
            #elif KB == 222
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 222
               movddup rB1, rB0
            #elif KB == 222
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 222
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 6992(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 6992(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 223
               movapd 7008(pB), rB1
            #elif KB == 223
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 223
               movapd 7008(pA), rA0
            #elif KB == 223
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 223
               movddup rB1, rB0
            #elif KB == 223
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 223
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7024(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7024(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 224
               movapd 7040(pB), rB1
            #elif KB == 224
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 224
               movapd 7040(pA), rA0
            #elif KB == 224
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 224
               movddup rB1, rB0
            #elif KB == 224
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 224
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7056(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7056(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 225
               movapd 7072(pB), rB1
            #elif KB == 225
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 225
               movapd 7072(pA), rA0
            #elif KB == 225
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 225
               movddup rB1, rB0
            #elif KB == 225
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 225
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7088(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7088(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 226
               movapd 7104(pB), rB1
            #elif KB == 226
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 226
               movapd 7104(pA), rA0
            #elif KB == 226
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 226
               movddup rB1, rB0
            #elif KB == 226
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 226
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7120(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7120(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 227
               movapd 7136(pB), rB1
            #elif KB == 227
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 227
               movapd 7136(pA), rA0
            #elif KB == 227
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 227
               movddup rB1, rB0
            #elif KB == 227
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 227
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7152(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7152(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 228
               movapd 7168(pB), rB1
            #elif KB == 228
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 228
               movapd 7168(pA), rA0
            #elif KB == 228
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 228
               movddup rB1, rB0
            #elif KB == 228
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 228
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7184(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7184(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 229
               movapd 7200(pB), rB1
            #elif KB == 229
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 229
               movapd 7200(pA), rA0
            #elif KB == 229
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 229
               movddup rB1, rB0
            #elif KB == 229
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 229
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7216(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7216(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 230
               movapd 7232(pB), rB1
            #elif KB == 230
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 230
               movapd 7232(pA), rA0
            #elif KB == 230
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 230
               movddup rB1, rB0
            #elif KB == 230
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 230
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7248(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7248(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 231
               movapd 7264(pB), rB1
            #elif KB == 231
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 231
               movapd 7264(pA), rA0
            #elif KB == 231
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 231
               movddup rB1, rB0
            #elif KB == 231
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 231
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7280(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7280(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 232
               movapd 7296(pB), rB1
            #elif KB == 232
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 232
               movapd 7296(pA), rA0
            #elif KB == 232
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 232
               movddup rB1, rB0
            #elif KB == 232
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 232
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7312(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7312(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 233
               movapd 7328(pB), rB1
            #elif KB == 233
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 233
               movapd 7328(pA), rA0
            #elif KB == 233
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 233
               movddup rB1, rB0
            #elif KB == 233
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 233
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7344(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7344(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 234
               movapd 7360(pB), rB1
            #elif KB == 234
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 234
               movapd 7360(pA), rA0
            #elif KB == 234
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 234
               movddup rB1, rB0
            #elif KB == 234
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 234
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7376(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7376(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 235
               movapd 7392(pB), rB1
            #elif KB == 235
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 235
               movapd 7392(pA), rA0
            #elif KB == 235
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 235
               movddup rB1, rB0
            #elif KB == 235
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 235
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7408(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7408(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 236
               movapd 7424(pB), rB1
            #elif KB == 236
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 236
               movapd 7424(pA), rA0
            #elif KB == 236
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 236
               movddup rB1, rB0
            #elif KB == 236
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 236
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7440(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7440(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 237
               movapd 7456(pB), rB1
            #elif KB == 237
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 237
               movapd 7456(pA), rA0
            #elif KB == 237
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 237
               movddup rB1, rB0
            #elif KB == 237
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 237
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7472(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7472(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 238
               movapd 7488(pB), rB1
            #elif KB == 238
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 238
               movapd 7488(pA), rA0
            #elif KB == 238
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 238
               movddup rB1, rB0
            #elif KB == 238
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 238
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7504(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7504(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 239
               movapd 7520(pB), rB1
            #elif KB == 239
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 239
               movapd 7520(pA), rA0
            #elif KB == 239
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 239
               movddup rB1, rB0
            #elif KB == 239
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 239
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7536(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7536(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 240
               movapd 7552(pB), rB1
            #elif KB == 240
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 240
               movapd 7552(pA), rA0
            #elif KB == 240
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 240
               movddup rB1, rB0
            #elif KB == 240
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 240
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7568(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7568(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 241
               movapd 7584(pB), rB1
            #elif KB == 241
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 241
               movapd 7584(pA), rA0
            #elif KB == 241
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 241
               movddup rB1, rB0
            #elif KB == 241
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 241
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7600(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7600(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 242
               movapd 7616(pB), rB1
            #elif KB == 242
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 242
               movapd 7616(pA), rA0
            #elif KB == 242
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 242
               movddup rB1, rB0
            #elif KB == 242
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 242
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7632(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7632(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 243
               movapd 7648(pB), rB1
            #elif KB == 243
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 243
               movapd 7648(pA), rA0
            #elif KB == 243
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 243
               movddup rB1, rB0
            #elif KB == 243
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 243
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7664(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7664(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 244
               movapd 7680(pB), rB1
            #elif KB == 244
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 244
               movapd 7680(pA), rA0
            #elif KB == 244
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 244
               movddup rB1, rB0
            #elif KB == 244
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 244
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7696(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7696(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 245
               movapd 7712(pB), rB1
            #elif KB == 245
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 245
               movapd 7712(pA), rA0
            #elif KB == 245
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 245
               movddup rB1, rB0
            #elif KB == 245
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 245
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7728(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7728(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 246
               movapd 7744(pB), rB1
            #elif KB == 246
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 246
               movapd 7744(pA), rA0
            #elif KB == 246
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 246
               movddup rB1, rB0
            #elif KB == 246
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 246
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7760(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7760(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 247
               movapd 7776(pB), rB1
            #elif KB == 247
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 247
               movapd 7776(pA), rA0
            #elif KB == 247
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 247
               movddup rB1, rB0
            #elif KB == 247
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 247
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7792(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7792(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 248
               movapd 7808(pB), rB1
            #elif KB == 248
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 248
               movapd 7808(pA), rA0
            #elif KB == 248
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 248
               movddup rB1, rB0
            #elif KB == 248
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 248
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7824(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7824(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 249
               movapd 7840(pB), rB1
            #elif KB == 249
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 249
               movapd 7840(pA), rA0
            #elif KB == 249
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 249
               movddup rB1, rB0
            #elif KB == 249
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 249
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7856(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7856(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 250
               movapd 7872(pB), rB1
            #elif KB == 250
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 250
               movapd 7872(pA), rA0
            #elif KB == 250
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 250
               movddup rB1, rB0
            #elif KB == 250
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 250
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7888(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7888(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 251
               movapd 7904(pB), rB1
            #elif KB == 251
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 251
               movapd 7904(pA), rA0
            #elif KB == 251
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 251
               movddup rB1, rB0
            #elif KB == 251
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 251
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7920(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7920(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 252
               movapd 7936(pB), rB1
            #elif KB == 252
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 252
               movapd 7936(pA), rA0
            #elif KB == 252
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 252
               movddup rB1, rB0
            #elif KB == 252
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 252
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7952(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7952(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 253
               movapd 7968(pB), rB1
            #elif KB == 253
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 253
               movapd 7968(pA), rA0
            #elif KB == 253
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 253
               movddup rB1, rB0
            #elif KB == 253
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 253
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 7984(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 7984(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 254
               movapd 8000(pB), rB1
            #elif KB == 254
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 254
               movapd 8000(pA), rA0
            #elif KB == 254
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 254
               movddup rB1, rB0
            #elif KB == 254
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 254
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 8016(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 8016(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 255
               movapd 8032(pB), rB1
            #elif KB == 255
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 255
               movapd 8032(pA), rA0
            #elif KB == 255
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 255
               movddup rB1, rB0
            #elif KB == 255
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if KB > 255
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 8048(pA), rA1
            mulpd rA1, rB0
            addpd rB0, rC10
            unpckhpd rB1, rB1

            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movapd 8048(pB), rB3
            mulpd rA1, rB1
            addpd rB1, rC11
            movddup rB3, rB2

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            unpckhpd rB3, rB3
            mulpd rA1, rB2
            addpd rB2, rC12

            #if KB > 256
               movapd 8064(pB), rB1
            #elif KB == 256
               pref (pf)
            #endif
            mulpd rB3, rA0
            addpd rA0, rC03
            #if KB > 256
               movapd 8064(pA), rA0
            #elif KB == 256
               pref (pfC)
               add incPF, pfC
            #endif
            mulpd rA1, rB3
            addpd rB3, rC13
            #if KB > 256
               movddup rB1, rB0
            #elif KB == 256
               prefetcht1 64(pfB)
               add incPF, pfB
               add incPF, pf
            #endif
         #endif
         #if defined(BETA1) || defined(BETAN1)
            BETCOP (pC), rC00
            movapd rC00, (pC)
            BETCOP 16(pC), rC10
            movapd rC10, 16(pC)
            BETCOP 32(pC), rC01
            movapd rC01, 32(pC)
            BETCOP 48(pC), rC11
            movapd rC11, 48(pC)
            BETCOP 64(pC), rC02
            movapd rC02, 64(pC)
            BETCOP 80(pC), rC12
            movapd rC12, 80(pC)
            BETCOP 96(pC), rC03
            movapd rC03, 96(pC)
            BETCOP 112(pC), rC13
            movapd rC13, 112(pC)
         #else
            movapd rC00, (pC)
            movapd rC10, 16(pC)
            movapd rC01, 32(pC)
            movapd rC11, 48(pC)
            movapd rC02, 64(pC)
            movapd rC12, 80(pC)
            movapd rC03, 96(pC)
            movapd rC13, 112(pC)
         #endif
         sub $-128, pC
         add $KB*4*8, pB
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
