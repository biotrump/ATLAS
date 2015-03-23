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
#define incB    %rax
#define KK      %rdx
/*
 * Second (N-) loop items get next level of priority on good regs
 */
#define pC      %rbp
#define pfA     %rsi
#define pfB     %rbx
#define K1      %r8
#define incPF   %r9
#define nnu     %r10
#define incAm   %r11
#define KK0     %r15
/*
 * Outer- (M-) loop variables assigned to any regs
 */
#define nmu     %r12
#define pB0     %r13
#define nnu0    %r14
/*
 * floating point registers
 */
#define a0      %xmm0
#define b0      %xmm1
#define rC0     %xmm2
#define rC1     %xmm3
#define rC2     %xmm4
#define rC3     %xmm5
#define rC4     %xmm6
#define rC5     %xmm7
#define rC6     %xmm8
#define rC7     %xmm9
#define rC8     %xmm10
#define rC9     %xmm11
#define rC10    %xmm12
#define rC11    %xmm13
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
   #define prefB prefetcht0
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
   mov $2*12*8, incPF           /* incPF = mu*nu*sizeof */
/*
 * Extend range of 1-byte offsets  by starting at -128
 */
   sub $-128, pB
   sub $-128, pC
   sub $-128, pfA
   sub $-128, pfB
   movq pB, pB0
#if KB == 0
/*
 * Make KK even, and store if we have an extra iteration in K1
 */
   mov $2, KK0  /* KK0 sets 2nd bit, used only when K==3 */
   mov $1, K1
   and KK, K1   /* K1 1 for odd K, else 0 */
   add K1, KK0  /* KK0 has both bits set appropriately for K==3 case */
   cmp $3, KK   /* is K == 3? */
   CMOVE KK0, K1/* bit 0: even/odd, bit1: set only if K == 3 */
   sub K1, KK   /* KK now guaranteed even */
#elif KB%2
      sub $1, KK  /* make loop variable even */
#endif
/*
 * incAm = KB*MU*sizeof = K*2*8 = K*16
 */
   shl $4, KK                   /* KK = K*MU*sizeof = K*2*8 = K*16 */
   mov KK, incAm                /* incAm = K*16 */
   mov $192, incB
   add KK, pA                   /* pA += K*MU (ptr starts 1 past last loc */
   neg KK
   mov KK, KK0

   nop ; nop ; nop ; nop ; nop
   .local MNLOOP
   MNLOOP:
/*
      .local NLOOP
      NLOOP:
*/
         #if KB == 0
            bt $1, K1   /* If K==3, use special-case code that avoids loop */
            jc KISTHREE
         #endif
/*
 *       Peel first iteration of K loop to initialize rCx
 */
         movapd (pA,KK), a0
         movddup -128(pB), rC0
         mulpd a0, rC0
            prefC -128(pC)
         movddup -120(pB), rC1
         mulpd a0, rC1
            prefC -64(pC)
         movddup -112(pB), rC2
         mulpd a0, rC2
            prefC (pC)
         movddup -104(pB), rC3
         mulpd a0, rC3
            prefA -128(pfA)
         movddup -96(pB), rC4
         mulpd a0, rC4
            prefA -64(pfA)
         movddup -88(pB), rC5
         mulpd a0, rC5
            prefA (pfA)
         movddup -80(pB), rC6
         mulpd a0, rC6
            prefB -128(pfB)
         movddup -72(pB), rC7
         mulpd a0, rC7
            prefB -64(pfB)
         movddup -64(pB), rC8
         mulpd a0, rC8
            prefB (pfB)
         movddup -56(pB), rC9
         mulpd a0, rC9
         movddup -48(pB), rC10
         mulpd a0, rC10
         movddup -40(pB), rC11
         mulpd a0, rC11
      #if KB == 0
         cmp $0, KK
         jz KLOOPDONE
      #endif

         #if KB > 1 || KB == 0
            movapd 16(pA,KK), a0
            movddup -32(pB), b0
            mulpd a0, b0
            addpd b0, rC0
            movddup -24(pB), b0
            mulpd a0, b0
            addpd b0, rC1
               #if KB > 2
                  add $32, KK
               #endif
            movddup -16(pB), b0
            mulpd a0, b0
            addpd b0, rC2
               add incPF, pfA
            movddup -8(pB), b0
            mulpd a0, b0
            addpd b0, rC3
               add incPF, pfB
            movddup (pB), b0
            mulpd a0, b0
            addpd b0, rC4
            movddup 8(pB), b0
            mulpd a0, b0
            addpd b0, rC5
            movddup 16(pB), b0
            mulpd a0, b0
            addpd b0, rC6
            movddup 24(pB), b0
            mulpd a0, b0
            addpd b0, rC7
            movddup 32(pB), b0
            mulpd a0, b0
            addpd b0, rC8
            movddup 40(pB), b0
            mulpd a0, b0
            addpd b0, rC9
            movddup 48(pB), b0
            mulpd a0, b0
            addpd b0, rC10
            movddup 56(pB), b0
            add incB, pB
            mulpd a0, b0
            addpd b0, rC11
            #if KB == 0
               add $32, KK
               jz KLOOPDONE
            #endif
         #endif

         #if KB >= 4 || KB == 0
         .local KLOOP
         KLOOP:
            movapd (pA,KK), a0
            movddup -128(pB), b0
            mulpd a0, b0
            addpd b0, rC0
            movddup -120(pB), b0
            mulpd a0, b0
            addpd b0, rC1
            movddup -112(pB), b0
            mulpd a0, b0
            addpd b0, rC2
            movddup -104(pB), b0
            mulpd a0, b0
            addpd b0, rC3
            movddup -96(pB), b0
            mulpd a0, b0
            addpd b0, rC4
            movddup -88(pB), b0
            mulpd a0, b0
            addpd b0, rC5
            movddup -80(pB), b0
            mulpd a0, b0
            addpd b0, rC6
            movddup -72(pB), b0
            mulpd a0, b0
            addpd b0, rC7
            movddup -64(pB), b0
            mulpd a0, b0
            addpd b0, rC8
            movddup -56(pB), b0
            mulpd a0, b0
            addpd b0, rC9
            movddup -48(pB), b0
            mulpd a0, b0
            addpd b0, rC10
            movddup -40(pB), b0
            mulpd a0, b0
            addpd b0, rC11

            movapd 16(pA,KK), a0
            movddup -32(pB), b0
            mulpd a0, b0
            addpd b0, rC0
            movddup -24(pB), b0
            mulpd a0, b0
            addpd b0, rC1
            movddup -16(pB), b0
            mulpd a0, b0
            addpd b0, rC2
            movddup -8(pB), b0
            mulpd a0, b0
            addpd b0, rC3
            movddup (pB), b0
            mulpd a0, b0
            addpd b0, rC4
            movddup 8(pB), b0
            mulpd a0, b0
            addpd b0, rC5
            movddup 16(pB), b0
            mulpd a0, b0
            addpd b0, rC6
            movddup 24(pB), b0
            mulpd a0, b0
            addpd b0, rC7
            movddup 32(pB), b0
            mulpd a0, b0
            addpd b0, rC8
            movddup 40(pB), b0
            mulpd a0, b0
            addpd b0, rC9
            movddup 48(pB), b0
            mulpd a0, b0
            addpd b0, rC10
            movddup 56(pB), b0
            mulpd a0, b0
            addpd b0, rC11
            add incB, pB
         add $32, KK
         jnz KLOOP
         #endif
         #if KB == 0
            bt $0, K1    /* set CF if odd # of K */
            jc ONEEXTRAKIT
            .local KLOOPDONE
            KLOOPDONE:
         #elif KB%2 && KB != 1
            movapd (pA,KK), a0
            movddup -128(pB), b0
            mulpd a0, b0
            addpd b0, rC0
            movddup -120(pB), b0
            mulpd a0, b0
            addpd b0, rC1
            movddup -112(pB), b0
            mulpd a0, b0
            addpd b0, rC2
            movddup -104(pB), b0
            mulpd a0, b0
            addpd b0, rC3
            movddup -96(pB), b0
            mulpd a0, b0
            addpd b0, rC4
            movddup -88(pB), b0
            mulpd a0, b0
            addpd b0, rC5
            movddup -80(pB), b0
            mulpd a0, b0
            addpd b0, rC6
            movddup -72(pB), b0
            mulpd a0, b0
            addpd b0, rC7
            movddup -64(pB), b0
            mulpd a0, b0
            addpd b0, rC8
            movddup -56(pB), b0
            mulpd a0, b0
            addpd b0, rC9
            movddup -48(pB), b0
            mulpd a0, b0
            addpd b0, rC10
            movddup -40(pB), b0
            mulpd a0, b0
            addpd b0, rC11
         #endif
/*
 *       Write answer back out to C
 */
         #ifdef BETA0
            movapd rC0, -128(pC)
            movapd rC1, -112(pC)
            movapd rC2, -96(pC)
            movapd rC3, -80(pC)
            movapd rC4, -64(pC)
            movapd rC5, -48(pC)
            movapd rC6, -32(pC)
            movapd rC7, -16(pC)
            movapd rC8, (pC)
            movapd rC9, 16(pC)
            movapd rC10, 32(pC)
            movapd rC11, 48(pC)
/*
 *          Add running sum in rCx with original C, then store back out
 */
         #else
            BETCOP -128(pC), rC0
            movapd rC0, -128(pC)
            BETCOP -112(pC), rC1
            movapd rC1, -112(pC)
            BETCOP -96(pC), rC2
            movapd rC2, -96(pC)
            BETCOP -80(pC), rC3
            movapd rC3, -80(pC)
            BETCOP -64(pC), rC4
            movapd rC4, -64(pC)
            BETCOP -48(pC), rC5
            movapd rC5, -48(pC)
            BETCOP -32(pC), rC6
            movapd rC6, -32(pC)
            BETCOP -16(pC), rC7
            movapd rC7, -16(pC)
            BETCOP (pC), rC8
            movapd rC8, (pC)
            BETCOP 16(pC), rC9
            movapd rC9, 16(pC)
            BETCOP 32(pC), rC10
            movapd rC10, 32(pC)
            BETCOP 48(pC), rC11
            movapd rC11, 48(pC)
         #endif
         mov KK0, KK
         add incB, pC   /* pC += 12*8*2 = 192 */
      sub $1, nnu
      jnz MNLOOP

      mov nnu0, nnu
      mov pB0, pB
      add incAm, pA    /* pA += KB*KU*MU*size */
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
/*
 * This code executed only for odd K
 */
#if KB == 0
KISTHREE:
/*
 *       Peel first 2 iteration of K to initialize rCx
 */
         movapd (pA,KK), a0
         movddup -128(pB), rC0
         mulpd a0, rC0
            prefC -128(pC)
         movddup -120(pB), rC1
         mulpd a0, rC1
            prefC -64(pC)
         movddup -112(pB), rC2
         mulpd a0, rC2
            prefC (pC)
         movddup -104(pB), rC3
         mulpd a0, rC3
            prefA -128(pfA)
         movddup -96(pB), rC4
         mulpd a0, rC4
            prefA -64(pfA)
         movddup -88(pB), rC5
         mulpd a0, rC5
            prefA (pfA)
         movddup -80(pB), rC6
         mulpd a0, rC6
            prefB -128(pfB)
         movddup -72(pB), rC7
         mulpd a0, rC7
            prefB -64(pfB)
         movddup -64(pB), rC8
         mulpd a0, rC8
            prefB (pfB)
         movddup -56(pB), rC9
         mulpd a0, rC9
         movddup -48(pB), rC10
         mulpd a0, rC10
         movddup -40(pB), rC11
         mulpd a0, rC11

         movapd 16(pA,KK), a0
         movddup -32(pB), b0
         mulpd a0, b0
         addpd b0, rC0
         movddup -24(pB), b0
         mulpd a0, b0
         addpd b0, rC1
         movddup -16(pB), b0
         mulpd a0, b0
         addpd b0, rC2
            add incPF, pfA
         movddup -8(pB), b0
         mulpd a0, b0
         addpd b0, rC3
            add incPF, pfB
         movddup (pB), b0
         mulpd a0, b0
         addpd b0, rC4
         movddup 8(pB), b0
         mulpd a0, b0
         addpd b0, rC5
         movddup 16(pB), b0
         mulpd a0, b0
         addpd b0, rC6
         movddup 24(pB), b0
         mulpd a0, b0
         addpd b0, rC7
         movddup 32(pB), b0
         mulpd a0, b0
         addpd b0, rC8
         movddup 40(pB), b0
         mulpd a0, b0
         addpd b0, rC9
         movddup 48(pB), b0
         add $32, KK
         mulpd a0, b0
         addpd b0, rC10
         movddup 56(pB), b0
         add incB, pB
         mulpd a0, b0
         addpd b0, rC11
/*
 *       Then fall thru into last K iteration to finish K==3
 */
.local ONEEXTRAKIT
ONEEXTRAKIT:
            movapd (pA,KK), a0
            movddup -128(pB), b0
            mulpd a0, b0
            addpd b0, rC0
            movddup -120(pB), b0
            mulpd a0, b0
            addpd b0, rC1
            movddup -112(pB), b0
            mulpd a0, b0
            addpd b0, rC2
            movddup -104(pB), b0
            mulpd a0, b0
            addpd b0, rC3
            movddup -96(pB), b0
            mulpd a0, b0
            addpd b0, rC4
            movddup -88(pB), b0
            mulpd a0, b0
            addpd b0, rC5
            movddup -80(pB), b0
            mulpd a0, b0
            addpd b0, rC6
            movddup -72(pB), b0
            mulpd a0, b0
            addpd b0, rC7
            movddup -64(pB), b0
            mulpd a0, b0
            addpd b0, rC8
            movddup -56(pB), b0
            mulpd a0, b0
            addpd b0, rC9
            movddup -48(pB), b0
            mulpd a0, b0
            addpd b0, rC10
            movddup -40(pB), b0
            mulpd a0, b0
            addpd b0, rC11
            jmp KLOOPDONE
#endif
