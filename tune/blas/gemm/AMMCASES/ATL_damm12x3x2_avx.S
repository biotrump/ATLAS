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
#if (KB/2)*2 != KB
   #error "KB must be a multiple of 2!"
#endif

/*
                    rdi      rsi    rdx        rcx         r8        r9
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                  8(%rsp)    16(%rsp)     24(%rsp)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */

#define FSIZE 6*8

#if defined(SREAL) || defined(SCPLX)
   #define vmovapd vmovaps
   #define vsubpd vsubps
   #define vaddpd vaddps
   #define vmulpd vmulps
   #define vbroadcastsd vbroadcastss
   #define SZ 4
#else
   #define vmovapd vmovaps
   #define SZ 8
#endif
#ifdef BETAN1
   #define VCOP vsubpd
#else
   #define VCOP vaddpd
#endif
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
/*
 * Then N-loop variables much less important, so use any orig regs left
 */
#define pC      %rbx  /* set after mov rsi to nnu () */
#define nnu     %rsi  /* comes in as rsi */
#define r128    %rbp
#define KK0     %rdi
/*
 * We could give a rat's ass about what registers used in outer (M-) loop
 */
#define nmu     %r8   /* comes in as rdi */
#define incAm   %r9
#define nnu0    %r10
#define pB0     %r11
/*
 * Prefetch definitions
 */
#define PFADIST 448
#define PFBDIST 128
#define prefA(pA_) prefetcht0 pA_
#define prefB(pB_) prefetcht0 pB_
#define prefC(pC_) prefetcht0 pC_
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
   nop
   nop
   movq %rbx, 8(%rsp)
#if 0
   movq %r12, 16(%rsp)
   movq %r13, 24(%rsp)
   movq %r14, 32(%rsp)
   movq %r15, 40(%rsp)
#endif
/*
 * Load paramaters
 */
   mov %r8, pB
   prefA(64(pA))
   mov %rdi, nmu
   mov $128, r128
   prefA(64(pA,r128))
   mov %r9, pC
   mov nnu, nnu0
   sub $-128, pC                 /* extend range of 1-byte offsets */
/*
 * Set constants
 */
   mov pB, pB0
/*
 * incAm = MU*sizeof*K = 12*8*K = 3*32*K
 */
   lea (KK, KK, 2), KK          /* KK = 3*K */
   shl $5, KK                   /* KK = 3*32*K = 12*8*K = MU*sizeof*K */
   mov KK, incAm                /* incAm = MU*sizeof*K */
   prefA((pA,r128,2))
   sub $3*2*32, KK              /* stop 1 iteration early for final peel */
   add KK, pA                   /* pA[-kk] will access A */
   neg KK                       /* KK = -NU*sizeof*K */
   mov KK, KK0
   vbroadcastsd (pB), rB0
   ALIGN32
   MNLOOP:
/*
 *          First peeled K=0 iteration gets to preloading next iter's data for
 *          loop while only doing 1 load/flop.  It does no adds to avoid
 *          having to zero the C registers.
 */

            vmovapd (pA,KK), rA0
            vmulpd rB0, rA0, rC00
            vmovapd 32(pA,KK), rA1
            vmulpd rB0, rA1, rC10
            vmovapd 64(pA,KK), rA2
            vmulpd rB0, rA2, rC20

            vbroadcastsd SZ(pB), rB1
            vmulpd rB1, rA0, rC01
            vbroadcastsd 2*SZ(pB), rB2
            vmulpd rB1, rA1, rC11
               vbroadcastsd 3*SZ(pB), rB0
            vmulpd rB1, rA2, rC21
               vbroadcastsd 4*SZ(pB), rB1

            vmulpd rB2, rA0, rC02
               vmovapd 96(pA,KK), rA0
            vmulpd rB2, rA1, rC12
               vmovapd 128(pA,KK), rA1
            vmulpd rB2, rA2, rC22
               vmovapd 160(pA,KK), rA2
               add $96, KK
/*
 *          For K=1 iteration of peel, prefetch C
 */
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5*SZ(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefC((pC))
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
               vbroadcastsd 6*SZ(pB), rB0

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
               prefC(64(pC))
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
               vbroadcastsd 7*SZ(pB), rB1
            add $6*SZ, pB
            add $96, KK

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
               vmovapd (pA,KK), rA0
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
               vmovapd 32(pA,KK), rA1
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
               vmovapd 64(pA,KK), rA2
/*
 *        Finally, start actual loop
 */
          ALIGN32
          KLOOP:
               prefB(PFBDIST(pB))
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2*SZ(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefA(PFADIST(pA,KK))
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
               vbroadcastsd 3*SZ(pB), rB0

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
               prefA(PFADIST+64(pA,KK))
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
               add $96, KK
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
               vbroadcastsd 4*SZ(pB), rB1

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
               vmovapd (pA,KK), rA0
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
               vmovapd 32(pA,KK), rA1
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
               vmovapd 64(pA,KK), rA2
/*
 *          K = 1 iteration
 */
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5*SZ(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefA(PFADIST+128(pA,KK))
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
               vbroadcastsd 6*SZ(pB), rB0

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
               prefB(3*KB*8(pB))
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
               vbroadcastsd 7*SZ(pB), rB1
            add $6*SZ, pB

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
               vmovapd 96(pA,KK), rA0
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
               vmovapd 96+32(pA,KK), rA1
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
               vmovapd 96+64(pA,KK), rA2
            add $96, KK
         jnz KLOOP
/*
 *       Last iteration peeled off bottom to allow store of C; this should
 *       strongly improve BETA=0, but may hurt BETA=1.  This kernel written
 *       primarily for K-cleanup, which is always BETA=0, so do it.
 */
         vmulpd rB0, rA0, m0
         vaddpd rC00, m0, rC00
         vbroadcastsd 2*SZ(pB), rB2
         vmulpd rB0, rA1, m0
         vaddpd rC10, m0, rC10
            prefA(PFADIST(pA))
         vmulpd rB0, rA2, m0
         vaddpd rC20, m0, rC20
            vbroadcastsd 3*SZ(pB), rB0

         vmulpd rB1, rA0, m0
         vaddpd rC01, m0, rC01
            prefA(PFADIST+64(pA))
         vmulpd rB1, rA1, m0
         vaddpd rC11, m0, rC11
         vmulpd rB1, rA2, m0
         vaddpd rC21, m0, rC21
            vbroadcastsd 4*SZ(pB), rB1

         vmulpd rB2, rA0, m0
         vaddpd rC02, m0, rC02
            vmovapd 96(pA), rA0
         vmulpd rB2, rA1, m0
         vaddpd rC12, m0, rC12
            vmovapd 96+32(pA), rA1
         vmulpd rB2, rA2, m0
         vaddpd rC22, m0, rC22
            vmovapd 96+64(pA), rA2
/*
 *       K = 1 iteration
 */
         vmulpd rB0, rA0, m0
         vaddpd rC00, m0, rC00
         #ifndef BETA0
            VCOP -128(pC), rC00, rC00
         #endif
         vmovapd rC00, -128(pC)
         vbroadcastsd 5*SZ(pB), rB2
         vmulpd rB0, rA1, m0
         vaddpd rC10, m0, rC10
         #ifndef BETA0
            VCOP -96(pC), rC10, rC10
         #endif
         vmovapd rC10, -96(pC)
            prefA(PFADIST+96+128(pA))
         vmulpd rB0, rA2, m0
         vaddpd rC20, m0, rC20
         #ifndef BETA0
            VCOP -64(pC), rC20, rC20
         #endif
         vmovapd rC20, -64(pC)
            vbroadcastsd 6*SZ(pB), rB0

         vmulpd rB1, rA0, m0
         vaddpd rC01, m0, rC01
         #ifndef BETA0
            VCOP -32(pC), rC01, rC01
         #endif
         vmovapd rC01, -32(pC)
            prefB(3*KB*8(pB))
         vmulpd rB1, rA1, m0
         vaddpd rC11, m0, rC11
         #ifndef BETA0
            VCOP (pC), rC11, rC11
         #endif
         vmovapd rC11, (pC)
         vmulpd rB1, rA2, m0
         vaddpd rC21, m0, rC21
         #ifndef BETA0
            VCOP 32(pC), rC21, rC21
         #endif
         vmovapd rC21, 32(pC)
         add $6*SZ, pB

         vmulpd rB2, rA0, m0
         vaddpd rC02, m0, rC02
         #ifndef BETA0
            VCOP 64(pC), rC02, rC02
         #endif
         vmovapd rC02, 64(pC)
         vmulpd rB2, rA1, m0
         vaddpd rC12, m0, rC12
         #ifndef BETA0
            VCOP 96(pC), rC12, rC12
         #endif
         vmovapd rC12, 96(pC)
         vmulpd rB2, rA2, m0
         vaddpd rC22, m0, rC22
         #ifndef BETA0
            VCOP (pC,r128), rC22, rC22
         #endif
         vmovapd rC22, (pC,r128)

         mov KK0, KK
         vbroadcastsd (pB), rB0
         add $12*3*8, pC                /* pC += MU*NU*sizeof */
      sub $1, nnu
      jnz MNLOOP
      vbroadcastsd (pB0), rB0
      mov nnu0, nnu
      mov pB0, pB
      add incAm, pA                     /* pA += MU*sizeof*K */
   sub $1, nmu
   jnz MNLOOP
#if 0
.local DONE
DONE:
#endif
   movq (%rsp), %rbp
   movq 8(%rsp), %rbx
#if 0
   movq 16(%rsp), %r12
   movq 24(%rsp), %r13
   movq 32(%rsp), %r14
   movq 40(%rsp), %r15
#endif
   add  $FSIZE, %rsp
   ret
