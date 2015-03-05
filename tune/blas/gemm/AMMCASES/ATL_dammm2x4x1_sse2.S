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
#define vmovapd vmovaps
#define nmu     %rdi
#define nnu     %rsi
#define nnu0    %r10
#define KK      %rdx
#define KK0     %r11
#define pA      %rcx
#define pB      %rax
#define pC      %r9
#define pf      %rbp
#define pB0     %r12
#define incPF   %rbx
#define pA0     %r13
#define incAm   %r14
#define pfB     %r15
#define FSIZE 6*8
#define FSIZE 6*8

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
   movq %r8, pB
   mov nnu, nnu0
   movq FSIZE+16(%rsp), pf      /* pf = pBn */
   movq FSIZE+8(%rsp), pfB      /* pfB = pAn */
   cmp pf, pB                   /* if (pBn == pB) */
   CMOVE pfB, pf                /* if (pBn == pB) pfB = pAn */
   CMOVEq FSIZE+24(%rsp), pfB   /* if (pbN == pB) pfB = pCn */
   mov $8*4*4, incPF
   mov KK, incAm                /* incAm = K */
   shl $5, incAm                /* incAm = K*MU*size = K*4*8 = K*32 = K << 5 */
   mov KK, KK0
   mov pA, pA0
   movq pB, pB0

   ALIGN8
   .local MNLOOP
   MNLOOP:
/*
 *       Peel first iteration of K-loop to handle init of C to 0
 */
         movapd (pB), rB1
         movddup rB1, rB0
         movapd (pA), rC00

         movapd rC00, rC01
         mulpd rB0, rC00
         unpckhpd rB1, rB1
         movapd 16(pA), rC10
         movapd rC10, rC11
         mulpd rB0, rC10

         movapd 16(pB), rB3
         movddup rB3, rB2
         movapd rC01, rC02
         mulpd rB1, rC01
         pref (pf)
         movapd rC11, rC12
         mulpd rB1, rC11
         unpckhpd rB3, rB3

         pref 64(pf)
         movapd rC02, rC03
         mulpd rB2, rC02
         movapd rC12, rC13
         mulpd rB2, rC12

         prefC (pC)
         mulpd rB3, rC03
         add $32, pB
         prefC 64(pC)
         mulpd rB3, rC13
         add incPF, pf
         add $32, pA
         sub $1, KK
         jz DONEK

         KLOOP:
            movapd (pB), rB1
            movddup rB1, rB0
            movapd (pA), rA0
            movapd rB0, rm0
            mulpd rA0, rm0
            addpd rm0, rC00
            movapd 16(pA), rA1
            movapd rB0, rm0
            mulpd rA1, rm0
            addpd rm0, rC10
            unpckhpd rB1, rB1

            movapd 16(pB), rB3
            movapd rB1, rm0
            mulpd rA0, rm0
            addpd rm0, rC01
            movddup rB3, rB2
            movapd rB1, rm0
            mulpd rA1, rm0
            addpd rm0, rC11
            unpckhpd rB3, rB3

            movapd rB2, rm0
            mulpd rA0, rm0
            addpd rm0, rC02
            add $32, pA
            movapd rB2, rm0
            mulpd rA1, rm0
            addpd rm0, rC12
            add $32, pB

            movapd rB3, rm0
            mulpd rA0, rm0
            addpd rm0, rC03
            movapd rB3, rm0
            mulpd rA1, rm0
            addpd rm0, rC13
         sub $1, KK
         jnz KLOOP
         DONEK:
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
         mov KK0, KK
         mov pA0, pA
      sub $1, nnu
      jnz MNLOOP
      mov nnu0, nnu
      mov pB0, pB
      add incAm, pA0
      mov pA0, pA
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
