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
   #define BETCOP vsubpd
#else
   #define BETCOP vaddpd
#endif
#define vmovapd vmovaps
/*
 * floating point registers
 */
#define m0   %ymm0
#define rA0  %ymm1
#define rA1  %ymm2
#define rA2  %ymm3
#define rB0  %ymm4
#define rB1  %ymm5
#define rB2  %ymm6
#define rC00 %ymm7
#define rC10 %ymm8
#define rC20 %ymm9
#define rC01 %ymm10
#define rC11 %ymm11
#define rC21 %ymm12
#define rC02 %ymm13
#define rC12 %ymm14
#define rC22 %ymm15
/*
 * Prioritize original registers for inner-loop operations, but inc regs
 * can be anything w/o changing opcode size, so use new regs for those
 */
#define KK      %rdx  /* API reg */
#define pA      %rcx  /* API reg */
#define pB      %rax  /* comes in as r9 */
#define incAk   %r9   /* set after mov r9 to pC () */
#define incBk   %r8   /* set after mov r8 to pB (rax) */
/*
 * Then N-loop variables much less important, so use any orig regs left
 */
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
.text
.global ATL_asmdecor(ATL_USERMM)
ALIGN16
ATL_asmdecor(ATL_USERMM):
/*
 * Save callee-saved iregs
 */
   sub  $FSIZE, %rsp
   movq %rbp, (%rsp)
   movq %rbx, 8(%rsp)
   movq %r12, 16(%rsp)
   movq %r13, 24(%rsp)
   movq %r14, 32(%rsp)
   movq %r15, 40(%rsp)
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
   sub $-128, pfA
   sub $-128, pfB
   sub $-128, pC                 /* extend range of 1-byte offsets */
/*
 * Set constants
 */
   mov $12*3*8, incPF           /* incPF = mu*nu*sizeof */
   mov $24, incBk               /* 24 = NU*sizeof = 3*8 = 24 */
   mov $96, incAk               /* 96 = MU*sizeof = 12*8 = 96 */
   mov pB, pB0
/*
 * incAm = MU*sizeof*K = 12*8*K = 3*32*K
 */
   lea (KK, KK, 2), KK          /* KK = 3*K */
   shl $5, KK                   /* KK = 3*32*K = 12*8*K = MU*sizeof*K */
   mov KK, incAm                /* incAm = MU*sizeof*K */
   add KK, pA                   /* pA[-kk] will access A */
   neg KK                       /* KK = -NU*sizeof*K */
   mov KK, KK0
   vbroadcastsd (pB), rB0
   MLOOP:
      NLOOP:
/*
 *          First peeled iteration gets us to preloading next iter's data for
 *          loop while only doing 1 load/flop.  It does no adds to avoid
 *          having to zero the C registers.
 */

            vmovapd (pA,KK), rA0
            vmulpd rB0, rA0, rC00
            vmovapd 32(pA,KK), rA1
            vmulpd rB0, rA1, rC10
            vmovapd 64(pA,KK), rA2
            vmulpd rB0, rA2, rC20

            vbroadcastsd 8(pB), rB1
            vmulpd rB1, rA0, rC01
            vbroadcastsd 16(pB), rB2
            vmulpd rB1, rA1, rC11
               vbroadcastsd 24(pB), rB0
            vmulpd rB1, rA2, rC21
               vbroadcastsd 32(pB), rB1

            vmulpd rB2, rA0, rC02
               vmovapd 96(pA,KK), rA0
            vmulpd rB2, rA1, rC12
               vmovapd 128(pA,KK), rA1
            vmulpd rB2, rA2, rC22
               vmovapd 160(pA,KK), rA2
            add $24, pB
            add $96, KK
            jz KDONE_NOFINALPEEL
/*
 *          Stop loop on iteration early to set up for C, so 2nd peel is
 *          from end (bottom) of loop
 */
            add $24, pB
            add $96, KK
            jz KLOOPDRAIN
/*
 *          Now peel a 3rd iteration (already peeled one from top and bottom)
 *          in order to do some prefetch.  This peel is the exact code from
 *          loop, with some prefetch commands added.
 */
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -8(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefA -128(pfA)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
               vbroadcastsd (pB), rB0

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
               prefA (pfA)
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
               vbroadcastsd 8(pB), rB1

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
               vmovapd (pA,KK), rA0
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
               vmovapd 32(pA,KK), rA1
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
               vmovapd 64(pA,KK), rA2
            add incPF, pfA
            add incPF, pfB
            add $24, pB
            add $96, KK
            jz KLOOPDRAIN
/*
 *        Finally, start actual loop
 */
          KLOOP:
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -8(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefetcht0 64(pB)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
               vbroadcastsd (pB), rB0

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
               prefetcht0 128(pA)
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
               vbroadcastsd 8(pB), rB1

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
               vmovapd (pA,KK), rA0
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
               vmovapd 32(pA,KK), rA1
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
               vmovapd 64(pA,KK), rA2
            add $24, pB
            add $96, KK
         jnz KLOOP
/*
 *       Last iteration peeled off bottom to allow store of C; this should
 *       strongly improve BETA=0, but may hurt BETA=1.  This kernel written
 *       primarily for K-cleanup, which is always BETA=0, so do it.
 */
.local KLOOPDRAIN
KLOOPDRAIN:
#ifdef BETAN1
   #define VCOP vsubpd
#else
   #define VCOP vaddpd
#endif
            vbroadcastsd -8(pB), rB2
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            #ifndef BETA0
               VCOP  -128(pC), rC00, rC00
            #endif
            vmovapd rC00, -128(pC)
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            #ifndef BETA0
               VCOP  -96(pC), rC10, rC10
            #endif
            vmovapd rC10, -96(pC)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #ifndef BETA0
               VCOP  -64(pC), rC20, rC20
            #endif
            vmovapd rC20, -64(pC)

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            #ifndef BETA0
               VCOP  -32(pC), rC01, rC01
            #endif
            vmovapd rC01, -32(pC)
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            #ifndef BETA0
               VCOP  (pC), rC11, rC11
            #endif
            vmovapd rC11, (pC)
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #ifndef BETA0
               VCOP  32(pC), rC21, rC21
            #endif
            vmovapd rC21, 32(pC)

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #ifndef BETA0
               VCOP  64(pC), rC02, rC02
            #endif
            vmovapd rC02, 64(pC)
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #ifndef BETA0
               VCOP  96(pC), rC12, rC12
            #endif
            vmovapd rC12, 96(pC)
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #ifndef BETA0
               VCOP  128(pC), rC22, rC22
            #endif
            vmovapd rC22, 128(pC)
.local KLOOPDONE
KLOOPDONE:
         mov KK0, KK
         vbroadcastsd (pB), rB0
         add $12*3*8, pC                /* pC += MU*NU*sizeof */
      sub $1, nnu
      jnz NLOOP
      vbroadcastsd (pB0), rB0
      mov nnu0, nnu
      mov pB0, pB
      add incAm, pA                     /* pA += MU*sizeof*K */
   sub $1, nmu
   jnz MLOOP
/* DONE: */
   movq (%rsp), %rbp
   movq 8(%rsp), %rbx
   movq 16(%rsp), %r12
   movq 24(%rsp), %r13
   movq 32(%rsp), %r14
   movq 40(%rsp), %r15
   add  $FSIZE, %rsp
   ret
/*
 * got answers in rCxx, just need to apply them to memory
 */
KDONE_NOFINALPEEL:
/*
 *       Write answer back out to C
 */
         #ifdef BETA0
            vmovapd rC00, -128(pC)
            vmovapd rC10, -96(pC)
            vmovapd rC20, -64(pC)
            vmovapd rC01, -32(pC)
            vmovapd rC11, (pC)
            vmovapd rC21, 32(pC)
            vmovapd rC02, 64(pC)
            vmovapd rC12, 96(pC)
            vmovapd rC22, 128(pC)
/*
 *          Add running sum in rCx with original C, then store back out
 */
         #else
            BETCOP -128(pC), rC00, rC00
            vmovapd rC00, -128(pC)
            BETCOP -96(pC), rC10, rC10
            vmovapd rC10, -96(pC)
            BETCOP -64(pC), rC20, rC20
            vmovapd rC20, -64(pC)
            BETCOP -32(pC), rC01, rC01
            vmovapd rC01, -32(pC)
            BETCOP (pC), rC11, rC11
            vmovapd rC11, (pC)
            BETCOP 32(pC), rC21, rC21
            vmovapd rC21, 32(pC)
            BETCOP 64(pC), rC02, rC02
            vmovapd rC02, 64(pC)
            BETCOP 96(pC), rC12, rC12
            vmovapd rC12, 96(pC)
            BETCOP 128(pC), rC22, rC22
            vmovapd rC22, 128(pC)
         #endif
         jmp KLOOPDONE
