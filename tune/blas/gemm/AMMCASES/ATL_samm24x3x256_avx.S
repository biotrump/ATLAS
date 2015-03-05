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
   #define BETCOP vsubps
#else
   #define BETCOP vaddps
#endif
#define vmovapd vmovaps
#define vmulpd vmulps
#define vaddpd vaddps
#define vbroadcastsd vbroadcastss
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
   mov $2*12*4, incPF           /* incPF = mu*nu*sizeof */
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
         vmovapd -128(pA), rA0
            vbroadcastss -128(pB), rB0
         vmulpd rA0, rB0, rC00
         vmovapd -96(pA), rA1
         vmulpd rA1, rB0, rC10
         vmovapd -64(pA), rA2
         vmulpd rA2, rB0, rC20
            vbroadcastss -124(pB), rB1
         vmulpd rA0, rB1, rC01
            vbroadcastss -120(pB), rB2
         vmulpd rA1, rB1, rC11
         #if KB > 1
            vbroadcastss -116(pB), rB0
         #endif
         vmulpd rA2, rB1, rC21
         #if KB > 1
            vbroadcastss -112(pB), rB1
         #endif
         vmulpd rA0, rB2, rC02
         #if KB > 1
            vmovapd -32(pA), rA0
         #endif
         vmulpd rA1, rB2, rC12
         #if KB > 1
            vmovapd (pA), rA1
         #endif
         vmulpd rA2, rB2, rC22
         #if KB > 1
            vmovapd 32(pA), rA2
         #endif
/*
 *       Fully unrolled K-loop
 */
         #if KB > 1
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -108(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefC (pC)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 2
            vbroadcastss -104(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
               prefC 64(pC)
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 2
            vbroadcastss -100(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 2
            vmovaps 64(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 2
            vmovaps 96(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 2
            vmovaps -128(pA,i256), rA2
            #endif

         #endif
         #if KB > 2
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -96(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefC 128(pC)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 3
            vbroadcastss -92(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 3
            vbroadcastss -88(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 3
            vmovaps -96(pA,i256), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 3
            vmovaps -64(pA,i256), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 3
            vmovaps -32(pA,i256), rA2
            #endif

         #endif
         #if KB > 3
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -84(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefA -128(pfA)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 4
            vbroadcastss -80(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 4
            vbroadcastss -76(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 4
            vmovaps (pA,i256), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 4
            vmovaps 32(pA,i256), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 4
            vmovaps 64(pA,i256), rA2
            #endif

         #endif
         #if KB > 4
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -72(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefA -64(pfA)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 5
            vbroadcastss -68(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 5
            vbroadcastss -64(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 5
            vmovaps 96(pA,i256), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 5
            vmovaps -128(pA,i256,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 5
            vmovaps -96(pA,i256,2), rA2
            #endif

         #endif
         #if KB > 5
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -60(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefA (pfA)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 6
            vbroadcastss -56(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 6
            vbroadcastss -52(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 6
            vmovaps -64(pA,i256,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 6
            vmovaps -32(pA,i256,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 6
            vmovaps (pA,i256,2), rA2
            #endif

         #endif
         #if KB > 6
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -48(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefB -128(pfB)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 7
            vbroadcastss -44(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 7
            vbroadcastss -40(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 7
            vmovaps 32(pA,i256,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 7
            vmovaps 64(pA,i256,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 7
            vmovaps 96(pA,i256,2), rA2
            #endif

         #endif
         #if KB > 7
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -36(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefB -64(pfB)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 8
            vbroadcastss -32(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 8
            vbroadcastss -28(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 8
            vmovaps -128(pA,i768), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 8
            vmovaps -96(pA,i768), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 8
            vmovaps -64(pA,i768), rA2
            #endif

         #endif
         #if KB > 8
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -24(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefB (pfB)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 9
            vbroadcastss -20(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 9
            vbroadcastss -16(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 9
            vmovaps -32(pA,i768), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 9
            vmovaps (pA,i768), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 9
            vmovaps 32(pA,i768), rA2
            #endif

         #endif
         #if KB > 9
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -12(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 10
            vbroadcastss -8(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 10
            vbroadcastss -4(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 10
            vmovaps 64(pA,i768), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 10
            vmovaps 96(pA,i768), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 10
            vmovaps -128(pA,i256,4), rA2
            #endif

         #endif
         #if KB > 10
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss (pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 11
            vbroadcastss 4(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 11
            vbroadcastss 8(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 11
            vmovaps -96(pA,i256,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 11
            vmovaps -64(pA,i256,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 11
            vmovaps -32(pA,i256,4), rA2
            #endif

         #endif
         #if KB > 11
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 12(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 12
            vbroadcastss 16(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 12
            vbroadcastss 20(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 12
            vmovaps (pA,i256,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 12
            vmovaps 32(pA,i256,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 12
            vmovaps 64(pA,i256,4), rA2
            #endif

         #endif
         #if KB > 12
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 24(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 13
            vbroadcastss 28(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 13
            vbroadcastss 32(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 13
            vmovaps 96(pA,i256,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 13
            vmovaps -128(pA,i1280), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 13
            vmovaps -96(pA,i1280), rA2
            #endif

         #endif
         #if KB > 13
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 36(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 14
            vbroadcastss 40(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 14
            vbroadcastss 44(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 14
            vmovaps -64(pA,i1280), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 14
            vmovaps -32(pA,i1280), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 14
            vmovaps (pA,i1280), rA2
            #endif

         #endif
         #if KB > 14
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 48(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 15
            vbroadcastss 52(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 15
            vbroadcastss 56(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 15
            vmovaps 32(pA,i1280), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 15
            vmovaps 64(pA,i1280), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 15
            vmovaps 96(pA,i1280), rA2
            #endif

         #endif
         #if KB > 15
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 60(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 16
            vbroadcastss 64(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 16
            vbroadcastss 68(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 16
            vmovaps -128(pA,i768,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 16
            vmovaps -96(pA,i768,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 16
            vmovaps -64(pA,i768,2), rA2
            #endif

         #endif
         #if KB > 16
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 72(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 17
            vbroadcastss 76(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 17
            vbroadcastss 80(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 17
            vmovaps -32(pA,i768,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 17
            vmovaps (pA,i768,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 17
            vmovaps 32(pA,i768,2), rA2
            #endif

         #endif
         #if KB > 17
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 84(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 18
            vbroadcastss 88(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 18
            vbroadcastss 92(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 18
            vmovaps 64(pA,i768,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 18
            vmovaps 96(pA,i768,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 18
            vmovaps -128(pA,i1792), rA2
            #endif

         #endif
         #if KB > 18
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 96(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 19
            vbroadcastss 100(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 19
            vbroadcastss 104(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 19
            vmovaps -96(pA,i1792), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 19
            vmovaps -64(pA,i1792), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 19
            vmovaps -32(pA,i1792), rA2
            #endif

         #endif
         #if KB > 19
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 108(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 20
            vbroadcastss 112(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 20
            vbroadcastss 116(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 20
            vmovaps (pA,i1792), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 20
            vmovaps 32(pA,i1792), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 20
            vmovaps 64(pA,i1792), rA2
            #endif

         #endif
         #if KB > 20
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 120(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 21
            vbroadcastss 124(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 21
            vbroadcastss -128(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 21
            vmovaps 96(pA,i1792), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 21
            vmovaps -128(pA,i256,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 21
            vmovaps -96(pA,i256,8), rA2
            #endif

         #endif
         #if KB > 21
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -124(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 22
            vbroadcastss -120(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 22
            vbroadcastss -116(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 22
            vmovaps -64(pA,i256,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 22
            vmovaps -32(pA,i256,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 22
            vmovaps (pA,i256,8), rA2
            #endif

         #endif
         #if KB > 22
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -112(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 23
            vbroadcastss -108(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 23
            vbroadcastss -104(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 23
            vmovaps 32(pA,i256,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 23
            vmovaps 64(pA,i256,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 23
            vmovaps 96(pA,i256,8), rA2
            #endif

         #endif
         #if KB > 23
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -100(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 24
            vbroadcastss -96(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 24
            vbroadcastss -92(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 24
            vmovaps -128(pA,i2304), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 24
            vmovaps -96(pA,i2304), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 24
            vmovaps -64(pA,i2304), rA2
            #endif

         #endif
         #if KB > 24
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -88(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 25
            vbroadcastss -84(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 25
            vbroadcastss -80(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 25
            vmovaps -32(pA,i2304), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 25
            vmovaps (pA,i2304), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 25
            vmovaps 32(pA,i2304), rA2
            #endif

         #endif
         #if KB > 25
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -76(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 26
            vbroadcastss -72(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 26
            vbroadcastss -68(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 26
            vmovaps 64(pA,i2304), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 26
            vmovaps 96(pA,i2304), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 26
            vmovaps -128(pA,i1280,2), rA2
            #endif

         #endif
         #if KB > 26
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -64(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 27
            vbroadcastss -60(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 27
            vbroadcastss -56(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 27
            vmovaps -96(pA,i1280,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 27
            vmovaps -64(pA,i1280,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 27
            vmovaps -32(pA,i1280,2), rA2
            #endif

         #endif
         #if KB > 27
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -52(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 28
            vbroadcastss -48(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 28
            vbroadcastss -44(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 28
            vmovaps (pA,i1280,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 28
            vmovaps 32(pA,i1280,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 28
            vmovaps 64(pA,i1280,2), rA2
            #endif

         #endif
         #if KB > 28
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -40(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 29
            vbroadcastss -36(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 29
            vbroadcastss -32(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 29
            vmovaps 96(pA,i1280,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 29
            vmovaps 2688(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 29
            vmovaps 2720(pA), rA2
            #endif

         #endif
         #if KB > 29
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -28(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 30
            vbroadcastss -24(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 30
            vbroadcastss -20(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 30
            vmovaps 2752(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 30
            vmovaps 2784(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 30
            vmovaps 2816(pA), rA2
            #endif

         #endif
         #if KB > 30
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -16(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 31
            vbroadcastss -12(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 31
            vbroadcastss -8(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 31
            vmovaps 2848(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 31
            vmovaps 2880(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 31
            vmovaps 2912(pA), rA2
            #endif

         #endif
         #if KB > 31
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -4(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 32
            vbroadcastss (pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 32
            vbroadcastss 4(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 32
            vmovaps -128(pA,i768,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 32
            vmovaps -96(pA,i768,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 32
            vmovaps -64(pA,i768,4), rA2
            #endif

         #endif
         #if KB > 32
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 8(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 33
            vbroadcastss 12(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 33
            vbroadcastss 16(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 33
            vmovaps -32(pA,i768,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 33
            vmovaps (pA,i768,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 33
            vmovaps 32(pA,i768,4), rA2
            #endif

         #endif
         #if KB > 33
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 20(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 34
            vbroadcastss 24(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 34
            vbroadcastss 28(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 34
            vmovaps 64(pA,i768,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 34
            vmovaps 96(pA,i768,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 34
            vmovaps 3200(pA), rA2
            #endif

         #endif
         #if KB > 34
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 32(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 35
            vbroadcastss 36(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 35
            vbroadcastss 40(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 35
            vmovaps 3232(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 35
            vmovaps 3264(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 35
            vmovaps 3296(pA), rA2
            #endif

         #endif
         #if KB > 35
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 44(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 36
            vbroadcastss 48(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 36
            vbroadcastss 52(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 36
            vmovaps 3328(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 36
            vmovaps 3360(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 36
            vmovaps 3392(pA), rA2
            #endif

         #endif
         #if KB > 36
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 56(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 37
            vbroadcastss 60(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 37
            vbroadcastss 64(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 37
            vmovaps 3424(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 37
            vmovaps -128(pA,i1792,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 37
            vmovaps -96(pA,i1792,2), rA2
            #endif

         #endif
         #if KB > 37
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 68(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 38
            vbroadcastss 72(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 38
            vbroadcastss 76(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 38
            vmovaps -64(pA,i1792,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 38
            vmovaps -32(pA,i1792,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 38
            vmovaps (pA,i1792,2), rA2
            #endif

         #endif
         #if KB > 38
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 80(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 39
            vbroadcastss 84(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 39
            vbroadcastss 88(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 39
            vmovaps 32(pA,i1792,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 39
            vmovaps 64(pA,i1792,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 39
            vmovaps 96(pA,i1792,2), rA2
            #endif

         #endif
         #if KB > 39
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 92(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 40
            vbroadcastss 96(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 40
            vbroadcastss 100(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 40
            vmovaps 3712(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 40
            vmovaps 3744(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 40
            vmovaps 3776(pA), rA2
            #endif

         #endif
         #if KB > 40
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 104(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 41
            vbroadcastss 108(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 41
            vbroadcastss 112(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 41
            vmovaps 3808(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 41
            vmovaps 3840(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 41
            vmovaps 3872(pA), rA2
            #endif

         #endif
         #if KB > 41
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 116(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 42
            vbroadcastss 120(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 42
            vbroadcastss 124(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 42
            vmovaps 3904(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 42
            vmovaps 3936(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 42
            vmovaps 3968(pA), rA2
            #endif

         #endif
         #if KB > 42
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -128(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 43
            vbroadcastss -124(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 43
            vbroadcastss -120(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 43
            vmovaps 4000(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 43
            vmovaps 4032(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 43
            vmovaps 4064(pA), rA2
            #endif

         #endif
         #if KB > 43
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -116(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 44
            vbroadcastss -112(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 44
            vbroadcastss -108(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 44
            vmovaps 4096(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 44
            vmovaps 4128(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 44
            vmovaps 4160(pA), rA2
            #endif

         #endif
         #if KB > 44
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -104(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 45
            vbroadcastss -100(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 45
            vbroadcastss -96(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 45
            vmovaps 4192(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 45
            vmovaps 4224(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 45
            vmovaps 4256(pA), rA2
            #endif

         #endif
         #if KB > 45
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -92(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 46
            vbroadcastss -88(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 46
            vbroadcastss -84(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 46
            vmovaps 4288(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 46
            vmovaps 4320(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 46
            vmovaps 4352(pA), rA2
            #endif

         #endif
         #if KB > 46
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -80(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 47
            vbroadcastss -76(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 47
            vbroadcastss -72(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 47
            vmovaps 4384(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 47
            vmovaps 4416(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 47
            vmovaps 4448(pA), rA2
            #endif

         #endif
         #if KB > 47
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -68(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 48
            vbroadcastss -64(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 48
            vbroadcastss -60(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 48
            vmovaps -128(pA,i2304,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 48
            vmovaps -96(pA,i2304,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 48
            vmovaps -64(pA,i2304,2), rA2
            #endif

         #endif
         #if KB > 48
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -56(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 49
            vbroadcastss -52(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 49
            vbroadcastss -48(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 49
            vmovaps -32(pA,i2304,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 49
            vmovaps (pA,i2304,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 49
            vmovaps 32(pA,i2304,2), rA2
            #endif

         #endif
         #if KB > 49
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -44(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 50
            vbroadcastss -40(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 50
            vbroadcastss -36(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 50
            vmovaps 64(pA,i2304,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 50
            vmovaps 96(pA,i2304,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 50
            vmovaps 4736(pA), rA2
            #endif

         #endif
         #if KB > 50
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -32(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 51
            vbroadcastss -28(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 51
            vbroadcastss -24(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 51
            vmovaps 4768(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 51
            vmovaps 4800(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 51
            vmovaps 4832(pA), rA2
            #endif

         #endif
         #if KB > 51
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -20(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 52
            vbroadcastss -16(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 52
            vbroadcastss -12(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 52
            vmovaps 4864(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 52
            vmovaps 4896(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 52
            vmovaps 4928(pA), rA2
            #endif

         #endif
         #if KB > 52
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -8(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 53
            vbroadcastss -4(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 53
            vbroadcastss (pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 53
            vmovaps 4960(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 53
            vmovaps -128(pA,i1280,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 53
            vmovaps -96(pA,i1280,4), rA2
            #endif

         #endif
         #if KB > 53
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 4(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 54
            vbroadcastss 8(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 54
            vbroadcastss 12(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 54
            vmovaps -64(pA,i1280,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 54
            vmovaps -32(pA,i1280,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 54
            vmovaps (pA,i1280,4), rA2
            #endif

         #endif
         #if KB > 54
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 16(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 55
            vbroadcastss 20(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 55
            vbroadcastss 24(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 55
            vmovaps 32(pA,i1280,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 55
            vmovaps 64(pA,i1280,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 55
            vmovaps 96(pA,i1280,4), rA2
            #endif

         #endif
         #if KB > 55
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 28(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 56
            vbroadcastss 32(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 56
            vbroadcastss 36(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 56
            vmovaps 5248(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 56
            vmovaps 5280(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 56
            vmovaps 5312(pA), rA2
            #endif

         #endif
         #if KB > 56
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 40(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 57
            vbroadcastss 44(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 57
            vbroadcastss 48(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 57
            vmovaps 5344(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 57
            vmovaps 5376(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 57
            vmovaps 5408(pA), rA2
            #endif

         #endif
         #if KB > 57
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 52(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 58
            vbroadcastss 56(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 58
            vbroadcastss 60(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 58
            vmovaps 5440(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 58
            vmovaps 5472(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 58
            vmovaps 5504(pA), rA2
            #endif

         #endif
         #if KB > 58
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 64(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 59
            vbroadcastss 68(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 59
            vbroadcastss 72(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 59
            vmovaps 5536(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 59
            vmovaps 5568(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 59
            vmovaps 5600(pA), rA2
            #endif

         #endif
         #if KB > 59
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 76(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 60
            vbroadcastss 80(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 60
            vbroadcastss 84(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 60
            vmovaps 5632(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 60
            vmovaps 5664(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 60
            vmovaps 5696(pA), rA2
            #endif

         #endif
         #if KB > 60
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 88(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 61
            vbroadcastss 92(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 61
            vbroadcastss 96(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 61
            vmovaps 5728(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 61
            vmovaps 5760(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 61
            vmovaps 5792(pA), rA2
            #endif

         #endif
         #if KB > 61
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 100(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 62
            vbroadcastss 104(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 62
            vbroadcastss 108(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 62
            vmovaps 5824(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 62
            vmovaps 5856(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 62
            vmovaps 5888(pA), rA2
            #endif

         #endif
         #if KB > 62
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 112(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 63
            vbroadcastss 116(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 63
            vbroadcastss 120(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 63
            vmovaps 5920(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 63
            vmovaps 5952(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 63
            vmovaps 5984(pA), rA2
            #endif

         #endif
         #if KB > 63
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 124(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 64
            vbroadcastss -128(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 64
            vbroadcastss -124(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 64
            vmovaps -128(pA,i768,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 64
            vmovaps -96(pA,i768,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 64
            vmovaps -64(pA,i768,8), rA2
            #endif

         #endif
         #if KB > 64
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -120(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 65
            vbroadcastss -116(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 65
            vbroadcastss -112(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 65
            vmovaps -32(pA,i768,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 65
            vmovaps (pA,i768,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 65
            vmovaps 32(pA,i768,8), rA2
            #endif

         #endif
         #if KB > 65
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -108(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 66
            vbroadcastss -104(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 66
            vbroadcastss -100(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 66
            vmovaps 64(pA,i768,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 66
            vmovaps 96(pA,i768,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 66
            vmovaps 6272(pA), rA2
            #endif

         #endif
         #if KB > 66
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -96(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 67
            vbroadcastss -92(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 67
            vbroadcastss -88(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 67
            vmovaps 6304(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 67
            vmovaps 6336(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 67
            vmovaps 6368(pA), rA2
            #endif

         #endif
         #if KB > 67
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -84(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 68
            vbroadcastss -80(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 68
            vbroadcastss -76(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 68
            vmovaps 6400(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 68
            vmovaps 6432(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 68
            vmovaps 6464(pA), rA2
            #endif

         #endif
         #if KB > 68
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -72(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 69
            vbroadcastss -68(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 69
            vbroadcastss -64(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 69
            vmovaps 6496(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 69
            vmovaps 6528(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 69
            vmovaps 6560(pA), rA2
            #endif

         #endif
         #if KB > 69
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -60(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 70
            vbroadcastss -56(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 70
            vbroadcastss -52(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 70
            vmovaps 6592(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 70
            vmovaps 6624(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 70
            vmovaps 6656(pA), rA2
            #endif

         #endif
         #if KB > 70
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -48(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 71
            vbroadcastss -44(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 71
            vbroadcastss -40(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 71
            vmovaps 6688(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 71
            vmovaps 6720(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 71
            vmovaps 6752(pA), rA2
            #endif

         #endif
         #if KB > 71
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -36(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 72
            vbroadcastss -32(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 72
            vbroadcastss -28(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 72
            vmovaps 6784(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 72
            vmovaps 6816(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 72
            vmovaps 6848(pA), rA2
            #endif

         #endif
         #if KB > 72
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -24(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 73
            vbroadcastss -20(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 73
            vbroadcastss -16(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 73
            vmovaps 6880(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 73
            vmovaps 6912(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 73
            vmovaps 6944(pA), rA2
            #endif

         #endif
         #if KB > 73
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -12(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 74
            vbroadcastss -8(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 74
            vbroadcastss -4(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 74
            vmovaps 6976(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 74
            vmovaps 7008(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 74
            vmovaps -128(pA,i1792,4), rA2
            #endif

         #endif
         #if KB > 74
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss (pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 75
            vbroadcastss 4(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 75
            vbroadcastss 8(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 75
            vmovaps -96(pA,i1792,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 75
            vmovaps -64(pA,i1792,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 75
            vmovaps -32(pA,i1792,4), rA2
            #endif

         #endif
         #if KB > 75
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 12(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 76
            vbroadcastss 16(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 76
            vbroadcastss 20(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 76
            vmovaps (pA,i1792,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 76
            vmovaps 32(pA,i1792,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 76
            vmovaps 64(pA,i1792,4), rA2
            #endif

         #endif
         #if KB > 76
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 24(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 77
            vbroadcastss 28(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 77
            vbroadcastss 32(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 77
            vmovaps 96(pA,i1792,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 77
            vmovaps 7296(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 77
            vmovaps 7328(pA), rA2
            #endif

         #endif
         #if KB > 77
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 36(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 78
            vbroadcastss 40(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 78
            vbroadcastss 44(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 78
            vmovaps 7360(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 78
            vmovaps 7392(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 78
            vmovaps 7424(pA), rA2
            #endif

         #endif
         #if KB > 78
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 48(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 79
            vbroadcastss 52(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 79
            vbroadcastss 56(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 79
            vmovaps 7456(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 79
            vmovaps 7488(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 79
            vmovaps 7520(pA), rA2
            #endif

         #endif
         #if KB > 79
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 60(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 80
            vbroadcastss 64(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 80
            vbroadcastss 68(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 80
            vmovaps 7552(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 80
            vmovaps 7584(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 80
            vmovaps 7616(pA), rA2
            #endif

         #endif
         #if KB > 80
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 72(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 81
            vbroadcastss 76(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 81
            vbroadcastss 80(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 81
            vmovaps 7648(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 81
            vmovaps 7680(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 81
            vmovaps 7712(pA), rA2
            #endif

         #endif
         #if KB > 81
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 84(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 82
            vbroadcastss 88(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 82
            vbroadcastss 92(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 82
            vmovaps 7744(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 82
            vmovaps 7776(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 82
            vmovaps 7808(pA), rA2
            #endif

         #endif
         #if KB > 82
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 96(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 83
            vbroadcastss 100(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 83
            vbroadcastss 104(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 83
            vmovaps 7840(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 83
            vmovaps 7872(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 83
            vmovaps 7904(pA), rA2
            #endif

         #endif
         #if KB > 83
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 108(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 84
            vbroadcastss 112(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 84
            vbroadcastss 116(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 84
            vmovaps 7936(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 84
            vmovaps 7968(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 84
            vmovaps 8000(pA), rA2
            #endif

         #endif
         #if KB > 84
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 120(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 85
            vbroadcastss 124(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 85
            vbroadcastss -128(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 85
            vmovaps 8032(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 85
            vmovaps 8064(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 85
            vmovaps 8096(pA), rA2
            #endif

         #endif
         #if KB > 85
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -124(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 86
            vbroadcastss -120(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 86
            vbroadcastss -116(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 86
            vmovaps 8128(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 86
            vmovaps 8160(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 86
            vmovaps 8192(pA), rA2
            #endif

         #endif
         #if KB > 86
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -112(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 87
            vbroadcastss -108(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 87
            vbroadcastss -104(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 87
            vmovaps 8224(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 87
            vmovaps 8256(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 87
            vmovaps 8288(pA), rA2
            #endif

         #endif
         #if KB > 87
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -100(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 88
            vbroadcastss -96(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 88
            vbroadcastss -92(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 88
            vmovaps 8320(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 88
            vmovaps 8352(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 88
            vmovaps 8384(pA), rA2
            #endif

         #endif
         #if KB > 88
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -88(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 89
            vbroadcastss -84(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 89
            vbroadcastss -80(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 89
            vmovaps 8416(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 89
            vmovaps 8448(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 89
            vmovaps 8480(pA), rA2
            #endif

         #endif
         #if KB > 89
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -76(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 90
            vbroadcastss -72(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 90
            vbroadcastss -68(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 90
            vmovaps 8512(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 90
            vmovaps 8544(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 90
            vmovaps 8576(pA), rA2
            #endif

         #endif
         #if KB > 90
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -64(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 91
            vbroadcastss -60(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 91
            vbroadcastss -56(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 91
            vmovaps 8608(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 91
            vmovaps 8640(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 91
            vmovaps 8672(pA), rA2
            #endif

         #endif
         #if KB > 91
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -52(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 92
            vbroadcastss -48(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 92
            vbroadcastss -44(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 92
            vmovaps 8704(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 92
            vmovaps 8736(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 92
            vmovaps 8768(pA), rA2
            #endif

         #endif
         #if KB > 92
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -40(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 93
            vbroadcastss -36(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 93
            vbroadcastss -32(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 93
            vmovaps 8800(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 93
            vmovaps 8832(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 93
            vmovaps 8864(pA), rA2
            #endif

         #endif
         #if KB > 93
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -28(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 94
            vbroadcastss -24(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 94
            vbroadcastss -20(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 94
            vmovaps 8896(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 94
            vmovaps 8928(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 94
            vmovaps 8960(pA), rA2
            #endif

         #endif
         #if KB > 94
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -16(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 95
            vbroadcastss -12(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 95
            vbroadcastss -8(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 95
            vmovaps 8992(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 95
            vmovaps 9024(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 95
            vmovaps 9056(pA), rA2
            #endif

         #endif
         #if KB > 95
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -4(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 96
            vbroadcastss (pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 96
            vbroadcastss 4(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 96
            vmovaps -128(pA,i2304,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 96
            vmovaps -96(pA,i2304,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 96
            vmovaps -64(pA,i2304,4), rA2
            #endif

         #endif
         #if KB > 96
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 8(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 97
            vbroadcastss 12(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 97
            vbroadcastss 16(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 97
            vmovaps -32(pA,i2304,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 97
            vmovaps (pA,i2304,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 97
            vmovaps 32(pA,i2304,4), rA2
            #endif

         #endif
         #if KB > 97
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 20(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 98
            vbroadcastss 24(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 98
            vbroadcastss 28(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 98
            vmovaps 64(pA,i2304,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 98
            vmovaps 96(pA,i2304,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 98
            vmovaps 9344(pA), rA2
            #endif

         #endif
         #if KB > 98
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 32(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 99
            vbroadcastss 36(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 99
            vbroadcastss 40(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 99
            vmovaps 9376(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 99
            vmovaps 9408(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 99
            vmovaps 9440(pA), rA2
            #endif

         #endif
         #if KB > 99
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 44(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 100
            vbroadcastss 48(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 100
            vbroadcastss 52(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 100
            vmovaps 9472(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 100
            vmovaps 9504(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 100
            vmovaps 9536(pA), rA2
            #endif

         #endif
         #if KB > 100
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 56(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 101
            vbroadcastss 60(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 101
            vbroadcastss 64(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 101
            vmovaps 9568(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 101
            vmovaps 9600(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 101
            vmovaps 9632(pA), rA2
            #endif

         #endif
         #if KB > 101
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 68(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 102
            vbroadcastss 72(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 102
            vbroadcastss 76(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 102
            vmovaps 9664(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 102
            vmovaps 9696(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 102
            vmovaps 9728(pA), rA2
            #endif

         #endif
         #if KB > 102
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 80(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 103
            vbroadcastss 84(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 103
            vbroadcastss 88(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 103
            vmovaps 9760(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 103
            vmovaps 9792(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 103
            vmovaps 9824(pA), rA2
            #endif

         #endif
         #if KB > 103
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 92(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 104
            vbroadcastss 96(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 104
            vbroadcastss 100(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 104
            vmovaps 9856(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 104
            vmovaps 9888(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 104
            vmovaps 9920(pA), rA2
            #endif

         #endif
         #if KB > 104
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 104(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 105
            vbroadcastss 108(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 105
            vbroadcastss 112(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 105
            vmovaps 9952(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 105
            vmovaps 9984(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 105
            vmovaps 10016(pA), rA2
            #endif

         #endif
         #if KB > 105
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 116(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 106
            vbroadcastss 120(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 106
            vbroadcastss 124(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 106
            vmovaps 10048(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 106
            vmovaps 10080(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 106
            vmovaps -128(pA,i1280,8), rA2
            #endif

         #endif
         #if KB > 106
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -128(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 107
            vbroadcastss -124(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 107
            vbroadcastss -120(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 107
            vmovaps -96(pA,i1280,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 107
            vmovaps -64(pA,i1280,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 107
            vmovaps -32(pA,i1280,8), rA2
            #endif

         #endif
         #if KB > 107
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -116(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 108
            vbroadcastss -112(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 108
            vbroadcastss -108(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 108
            vmovaps (pA,i1280,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 108
            vmovaps 32(pA,i1280,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 108
            vmovaps 64(pA,i1280,8), rA2
            #endif

         #endif
         #if KB > 108
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -104(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 109
            vbroadcastss -100(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 109
            vbroadcastss -96(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 109
            vmovaps 96(pA,i1280,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 109
            vmovaps 10368(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 109
            vmovaps 10400(pA), rA2
            #endif

         #endif
         #if KB > 109
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -92(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 110
            vbroadcastss -88(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 110
            vbroadcastss -84(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 110
            vmovaps 10432(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 110
            vmovaps 10464(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 110
            vmovaps 10496(pA), rA2
            #endif

         #endif
         #if KB > 110
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -80(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 111
            vbroadcastss -76(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 111
            vbroadcastss -72(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 111
            vmovaps 10528(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 111
            vmovaps 10560(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 111
            vmovaps 10592(pA), rA2
            #endif

         #endif
         #if KB > 111
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -68(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 112
            vbroadcastss -64(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 112
            vbroadcastss -60(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 112
            vmovaps 10624(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 112
            vmovaps 10656(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 112
            vmovaps 10688(pA), rA2
            #endif

         #endif
         #if KB > 112
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -56(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 113
            vbroadcastss -52(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 113
            vbroadcastss -48(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 113
            vmovaps 10720(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 113
            vmovaps 10752(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 113
            vmovaps 10784(pA), rA2
            #endif

         #endif
         #if KB > 113
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -44(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 114
            vbroadcastss -40(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 114
            vbroadcastss -36(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 114
            vmovaps 10816(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 114
            vmovaps 10848(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 114
            vmovaps 10880(pA), rA2
            #endif

         #endif
         #if KB > 114
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -32(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 115
            vbroadcastss -28(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 115
            vbroadcastss -24(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 115
            vmovaps 10912(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 115
            vmovaps 10944(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 115
            vmovaps 10976(pA), rA2
            #endif

         #endif
         #if KB > 115
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -20(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 116
            vbroadcastss -16(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 116
            vbroadcastss -12(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 116
            vmovaps 11008(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 116
            vmovaps 11040(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 116
            vmovaps 11072(pA), rA2
            #endif

         #endif
         #if KB > 116
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -8(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 117
            vbroadcastss -4(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 117
            vbroadcastss (pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 117
            vmovaps 11104(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 117
            vmovaps 11136(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 117
            vmovaps 11168(pA), rA2
            #endif

         #endif
         #if KB > 117
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 4(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 118
            vbroadcastss 8(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 118
            vbroadcastss 12(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 118
            vmovaps 11200(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 118
            vmovaps 11232(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 118
            vmovaps 11264(pA), rA2
            #endif

         #endif
         #if KB > 118
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 16(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 119
            vbroadcastss 20(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 119
            vbroadcastss 24(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 119
            vmovaps 11296(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 119
            vmovaps 11328(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 119
            vmovaps 11360(pA), rA2
            #endif

         #endif
         #if KB > 119
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 28(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 120
            vbroadcastss 32(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 120
            vbroadcastss 36(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 120
            vmovaps 11392(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 120
            vmovaps 11424(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 120
            vmovaps 11456(pA), rA2
            #endif

         #endif
         #if KB > 120
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 40(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 121
            vbroadcastss 44(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 121
            vbroadcastss 48(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 121
            vmovaps 11488(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 121
            vmovaps 11520(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 121
            vmovaps 11552(pA), rA2
            #endif

         #endif
         #if KB > 121
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 52(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 122
            vbroadcastss 56(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 122
            vbroadcastss 60(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 122
            vmovaps 11584(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 122
            vmovaps 11616(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 122
            vmovaps 11648(pA), rA2
            #endif

         #endif
         #if KB > 122
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 64(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 123
            vbroadcastss 68(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 123
            vbroadcastss 72(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 123
            vmovaps 11680(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 123
            vmovaps 11712(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 123
            vmovaps 11744(pA), rA2
            #endif

         #endif
         #if KB > 123
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 76(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 124
            vbroadcastss 80(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 124
            vbroadcastss 84(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 124
            vmovaps 11776(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 124
            vmovaps 11808(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 124
            vmovaps 11840(pA), rA2
            #endif

         #endif
         #if KB > 124
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 88(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 125
            vbroadcastss 92(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 125
            vbroadcastss 96(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 125
            vmovaps 11872(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 125
            vmovaps 11904(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 125
            vmovaps 11936(pA), rA2
            #endif

         #endif
         #if KB > 125
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 100(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 126
            vbroadcastss 104(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 126
            vbroadcastss 108(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 126
            vmovaps 11968(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 126
            vmovaps 12000(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 126
            vmovaps 12032(pA), rA2
            #endif

         #endif
         #if KB > 126
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 112(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 127
            vbroadcastss 116(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 127
            vbroadcastss 120(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 127
            vmovaps 12064(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 127
            vmovaps 12096(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 127
            vmovaps 12128(pA), rA2
            #endif

         #endif
         #if KB > 127
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 124(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 128
            vbroadcastss -128(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 128
            vbroadcastss -124(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 128
            vmovaps 12160(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 128
            vmovaps 12192(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 128
            vmovaps 12224(pA), rA2
            #endif

         #endif
         #if KB > 128
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -120(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 129
            vbroadcastss -116(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 129
            vbroadcastss -112(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 129
            vmovaps 12256(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 129
            vmovaps 12288(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 129
            vmovaps 12320(pA), rA2
            #endif

         #endif
         #if KB > 129
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -108(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 130
            vbroadcastss -104(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 130
            vbroadcastss -100(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 130
            vmovaps 12352(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 130
            vmovaps 12384(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 130
            vmovaps 12416(pA), rA2
            #endif

         #endif
         #if KB > 130
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -96(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 131
            vbroadcastss -92(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 131
            vbroadcastss -88(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 131
            vmovaps 12448(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 131
            vmovaps 12480(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 131
            vmovaps 12512(pA), rA2
            #endif

         #endif
         #if KB > 131
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -84(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 132
            vbroadcastss -80(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 132
            vbroadcastss -76(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 132
            vmovaps 12544(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 132
            vmovaps 12576(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 132
            vmovaps 12608(pA), rA2
            #endif

         #endif
         #if KB > 132
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -72(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 133
            vbroadcastss -68(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 133
            vbroadcastss -64(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 133
            vmovaps 12640(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 133
            vmovaps 12672(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 133
            vmovaps 12704(pA), rA2
            #endif

         #endif
         #if KB > 133
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -60(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 134
            vbroadcastss -56(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 134
            vbroadcastss -52(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 134
            vmovaps 12736(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 134
            vmovaps 12768(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 134
            vmovaps 12800(pA), rA2
            #endif

         #endif
         #if KB > 134
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -48(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 135
            vbroadcastss -44(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 135
            vbroadcastss -40(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 135
            vmovaps 12832(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 135
            vmovaps 12864(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 135
            vmovaps 12896(pA), rA2
            #endif

         #endif
         #if KB > 135
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -36(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 136
            vbroadcastss -32(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 136
            vbroadcastss -28(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 136
            vmovaps 12928(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 136
            vmovaps 12960(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 136
            vmovaps 12992(pA), rA2
            #endif

         #endif
         #if KB > 136
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -24(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 137
            vbroadcastss -20(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 137
            vbroadcastss -16(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 137
            vmovaps 13024(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 137
            vmovaps 13056(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 137
            vmovaps 13088(pA), rA2
            #endif

         #endif
         #if KB > 137
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -12(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 138
            vbroadcastss -8(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 138
            vbroadcastss -4(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 138
            vmovaps 13120(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 138
            vmovaps 13152(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 138
            vmovaps 13184(pA), rA2
            #endif

         #endif
         #if KB > 138
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss (pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 139
            vbroadcastss 4(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 139
            vbroadcastss 8(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 139
            vmovaps 13216(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 139
            vmovaps 13248(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 139
            vmovaps 13280(pA), rA2
            #endif

         #endif
         #if KB > 139
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 12(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 140
            vbroadcastss 16(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 140
            vbroadcastss 20(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 140
            vmovaps 13312(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 140
            vmovaps 13344(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 140
            vmovaps 13376(pA), rA2
            #endif

         #endif
         #if KB > 140
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 24(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 141
            vbroadcastss 28(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 141
            vbroadcastss 32(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 141
            vmovaps 13408(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 141
            vmovaps 13440(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 141
            vmovaps 13472(pA), rA2
            #endif

         #endif
         #if KB > 141
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 36(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 142
            vbroadcastss 40(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 142
            vbroadcastss 44(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 142
            vmovaps 13504(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 142
            vmovaps 13536(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 142
            vmovaps 13568(pA), rA2
            #endif

         #endif
         #if KB > 142
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 48(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 143
            vbroadcastss 52(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 143
            vbroadcastss 56(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 143
            vmovaps 13600(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 143
            vmovaps 13632(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 143
            vmovaps 13664(pA), rA2
            #endif

         #endif
         #if KB > 143
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 60(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 144
            vbroadcastss 64(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 144
            vbroadcastss 68(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 144
            vmovaps 13696(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 144
            vmovaps 13728(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 144
            vmovaps 13760(pA), rA2
            #endif

         #endif
         #if KB > 144
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 72(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 145
            vbroadcastss 76(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 145
            vbroadcastss 80(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 145
            vmovaps 13792(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 145
            vmovaps 13824(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 145
            vmovaps 13856(pA), rA2
            #endif

         #endif
         #if KB > 145
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 84(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 146
            vbroadcastss 88(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 146
            vbroadcastss 92(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 146
            vmovaps 13888(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 146
            vmovaps 13920(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 146
            vmovaps 13952(pA), rA2
            #endif

         #endif
         #if KB > 146
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 96(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 147
            vbroadcastss 100(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 147
            vbroadcastss 104(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 147
            vmovaps 13984(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 147
            vmovaps 14016(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 147
            vmovaps 14048(pA), rA2
            #endif

         #endif
         #if KB > 147
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 108(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 148
            vbroadcastss 112(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 148
            vbroadcastss 116(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 148
            vmovaps 14080(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 148
            vmovaps 14112(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 148
            vmovaps 14144(pA), rA2
            #endif

         #endif
         #if KB > 148
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 120(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 149
            vbroadcastss 124(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 149
            vbroadcastss -128(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 149
            vmovaps 14176(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 149
            vmovaps 14208(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 149
            vmovaps 14240(pA), rA2
            #endif

         #endif
         #if KB > 149
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -124(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 150
            vbroadcastss -120(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 150
            vbroadcastss -116(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 150
            vmovaps 14272(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 150
            vmovaps 14304(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 150
            vmovaps 14336(pA), rA2
            #endif

         #endif
         #if KB > 150
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -112(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 151
            vbroadcastss -108(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 151
            vbroadcastss -104(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 151
            vmovaps 14368(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 151
            vmovaps 14400(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 151
            vmovaps 14432(pA), rA2
            #endif

         #endif
         #if KB > 151
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -100(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 152
            vbroadcastss -96(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 152
            vbroadcastss -92(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 152
            vmovaps 14464(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 152
            vmovaps 14496(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 152
            vmovaps 14528(pA), rA2
            #endif

         #endif
         #if KB > 152
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -88(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 153
            vbroadcastss -84(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 153
            vbroadcastss -80(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 153
            vmovaps 14560(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 153
            vmovaps 14592(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 153
            vmovaps 14624(pA), rA2
            #endif

         #endif
         #if KB > 153
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -76(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 154
            vbroadcastss -72(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 154
            vbroadcastss -68(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 154
            vmovaps 14656(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 154
            vmovaps 14688(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 154
            vmovaps 14720(pA), rA2
            #endif

         #endif
         #if KB > 154
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -64(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 155
            vbroadcastss -60(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 155
            vbroadcastss -56(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 155
            vmovaps 14752(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 155
            vmovaps 14784(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 155
            vmovaps 14816(pA), rA2
            #endif

         #endif
         #if KB > 155
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -52(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 156
            vbroadcastss -48(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 156
            vbroadcastss -44(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 156
            vmovaps 14848(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 156
            vmovaps 14880(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 156
            vmovaps 14912(pA), rA2
            #endif

         #endif
         #if KB > 156
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -40(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 157
            vbroadcastss -36(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 157
            vbroadcastss -32(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 157
            vmovaps 14944(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 157
            vmovaps 14976(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 157
            vmovaps 15008(pA), rA2
            #endif

         #endif
         #if KB > 157
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -28(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 158
            vbroadcastss -24(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 158
            vbroadcastss -20(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 158
            vmovaps 15040(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 158
            vmovaps 15072(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 158
            vmovaps 15104(pA), rA2
            #endif

         #endif
         #if KB > 158
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -16(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 159
            vbroadcastss -12(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 159
            vbroadcastss -8(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 159
            vmovaps 15136(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 159
            vmovaps 15168(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 159
            vmovaps 15200(pA), rA2
            #endif

         #endif
         #if KB > 159
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -4(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 160
            vbroadcastss (pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 160
            vbroadcastss 4(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 160
            vmovaps 15232(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 160
            vmovaps 15264(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 160
            vmovaps 15296(pA), rA2
            #endif

         #endif
         #if KB > 160
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 8(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 161
            vbroadcastss 12(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 161
            vbroadcastss 16(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 161
            vmovaps 15328(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 161
            vmovaps 15360(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 161
            vmovaps 15392(pA), rA2
            #endif

         #endif
         #if KB > 161
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 20(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 162
            vbroadcastss 24(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 162
            vbroadcastss 28(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 162
            vmovaps 15424(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 162
            vmovaps 15456(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 162
            vmovaps 15488(pA), rA2
            #endif

         #endif
         #if KB > 162
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 32(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 163
            vbroadcastss 36(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 163
            vbroadcastss 40(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 163
            vmovaps 15520(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 163
            vmovaps 15552(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 163
            vmovaps 15584(pA), rA2
            #endif

         #endif
         #if KB > 163
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 44(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 164
            vbroadcastss 48(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 164
            vbroadcastss 52(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 164
            vmovaps 15616(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 164
            vmovaps 15648(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 164
            vmovaps 15680(pA), rA2
            #endif

         #endif
         #if KB > 164
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 56(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 165
            vbroadcastss 60(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 165
            vbroadcastss 64(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 165
            vmovaps 15712(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 165
            vmovaps 15744(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 165
            vmovaps 15776(pA), rA2
            #endif

         #endif
         #if KB > 165
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 68(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 166
            vbroadcastss 72(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 166
            vbroadcastss 76(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 166
            vmovaps 15808(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 166
            vmovaps 15840(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 166
            vmovaps 15872(pA), rA2
            #endif

         #endif
         #if KB > 166
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 80(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 167
            vbroadcastss 84(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 167
            vbroadcastss 88(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 167
            vmovaps 15904(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 167
            vmovaps 15936(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 167
            vmovaps 15968(pA), rA2
            #endif

         #endif
         #if KB > 167
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 92(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 168
            vbroadcastss 96(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 168
            vbroadcastss 100(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 168
            vmovaps 16000(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 168
            vmovaps 16032(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 168
            vmovaps 16064(pA), rA2
            #endif

         #endif
         #if KB > 168
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 104(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 169
            vbroadcastss 108(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 169
            vbroadcastss 112(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 169
            vmovaps 16096(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 169
            vmovaps 16128(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 169
            vmovaps 16160(pA), rA2
            #endif

         #endif
         #if KB > 169
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 116(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 170
            vbroadcastss 120(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 170
            vbroadcastss 124(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 170
            vmovaps 16192(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 170
            vmovaps 16224(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 170
            vmovaps 16256(pA), rA2
            #endif

         #endif
         #if KB > 170
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -128(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 171
            vbroadcastss -124(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 171
            vbroadcastss -120(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 171
            vmovaps 16288(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 171
            vmovaps 16320(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 171
            vmovaps 16352(pA), rA2
            #endif

         #endif
         #if KB > 171
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -116(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 172
            vbroadcastss -112(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 172
            vbroadcastss -108(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 172
            vmovaps 16384(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 172
            vmovaps 16416(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 172
            vmovaps 16448(pA), rA2
            #endif

         #endif
         #if KB > 172
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -104(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 173
            vbroadcastss -100(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 173
            vbroadcastss -96(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 173
            vmovaps 16480(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 173
            vmovaps 16512(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 173
            vmovaps 16544(pA), rA2
            #endif

         #endif
         #if KB > 173
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -92(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 174
            vbroadcastss -88(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 174
            vbroadcastss -84(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 174
            vmovaps 16576(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 174
            vmovaps 16608(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 174
            vmovaps 16640(pA), rA2
            #endif

         #endif
         #if KB > 174
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -80(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 175
            vbroadcastss -76(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 175
            vbroadcastss -72(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 175
            vmovaps 16672(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 175
            vmovaps 16704(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 175
            vmovaps 16736(pA), rA2
            #endif

         #endif
         #if KB > 175
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -68(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 176
            vbroadcastss -64(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 176
            vbroadcastss -60(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 176
            vmovaps 16768(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 176
            vmovaps 16800(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 176
            vmovaps 16832(pA), rA2
            #endif

         #endif
         #if KB > 176
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -56(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 177
            vbroadcastss -52(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 177
            vbroadcastss -48(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 177
            vmovaps 16864(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 177
            vmovaps 16896(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 177
            vmovaps 16928(pA), rA2
            #endif

         #endif
         #if KB > 177
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -44(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 178
            vbroadcastss -40(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 178
            vbroadcastss -36(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 178
            vmovaps 16960(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 178
            vmovaps 16992(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 178
            vmovaps 17024(pA), rA2
            #endif

         #endif
         #if KB > 178
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -32(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 179
            vbroadcastss -28(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 179
            vbroadcastss -24(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 179
            vmovaps 17056(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 179
            vmovaps 17088(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 179
            vmovaps 17120(pA), rA2
            #endif

         #endif
         #if KB > 179
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -20(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 180
            vbroadcastss -16(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 180
            vbroadcastss -12(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 180
            vmovaps 17152(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 180
            vmovaps 17184(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 180
            vmovaps 17216(pA), rA2
            #endif

         #endif
         #if KB > 180
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -8(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 181
            vbroadcastss -4(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 181
            vbroadcastss (pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 181
            vmovaps 17248(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 181
            vmovaps 17280(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 181
            vmovaps 17312(pA), rA2
            #endif

         #endif
         #if KB > 181
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 4(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 182
            vbroadcastss 8(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 182
            vbroadcastss 12(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 182
            vmovaps 17344(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 182
            vmovaps 17376(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 182
            vmovaps 17408(pA), rA2
            #endif

         #endif
         #if KB > 182
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 16(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 183
            vbroadcastss 20(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 183
            vbroadcastss 24(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 183
            vmovaps 17440(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 183
            vmovaps 17472(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 183
            vmovaps 17504(pA), rA2
            #endif

         #endif
         #if KB > 183
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 28(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 184
            vbroadcastss 32(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 184
            vbroadcastss 36(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 184
            vmovaps 17536(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 184
            vmovaps 17568(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 184
            vmovaps 17600(pA), rA2
            #endif

         #endif
         #if KB > 184
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 40(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 185
            vbroadcastss 44(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 185
            vbroadcastss 48(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 185
            vmovaps 17632(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 185
            vmovaps 17664(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 185
            vmovaps 17696(pA), rA2
            #endif

         #endif
         #if KB > 185
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 52(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 186
            vbroadcastss 56(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 186
            vbroadcastss 60(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 186
            vmovaps 17728(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 186
            vmovaps 17760(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 186
            vmovaps 17792(pA), rA2
            #endif

         #endif
         #if KB > 186
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 64(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 187
            vbroadcastss 68(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 187
            vbroadcastss 72(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 187
            vmovaps 17824(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 187
            vmovaps 17856(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 187
            vmovaps 17888(pA), rA2
            #endif

         #endif
         #if KB > 187
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 76(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 188
            vbroadcastss 80(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 188
            vbroadcastss 84(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 188
            vmovaps 17920(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 188
            vmovaps 17952(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 188
            vmovaps 17984(pA), rA2
            #endif

         #endif
         #if KB > 188
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 88(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 189
            vbroadcastss 92(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 189
            vbroadcastss 96(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 189
            vmovaps 18016(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 189
            vmovaps 18048(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 189
            vmovaps 18080(pA), rA2
            #endif

         #endif
         #if KB > 189
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 100(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 190
            vbroadcastss 104(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 190
            vbroadcastss 108(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 190
            vmovaps 18112(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 190
            vmovaps 18144(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 190
            vmovaps 18176(pA), rA2
            #endif

         #endif
         #if KB > 190
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 112(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 191
            vbroadcastss 116(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 191
            vbroadcastss 120(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 191
            vmovaps 18208(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 191
            vmovaps 18240(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 191
            vmovaps 18272(pA), rA2
            #endif

         #endif
         #if KB > 191
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 124(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 192
            vbroadcastss -128(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 192
            vbroadcastss -124(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 192
            vmovaps 18304(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 192
            vmovaps 18336(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 192
            vmovaps 18368(pA), rA2
            #endif

         #endif
         #if KB > 192
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -120(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 193
            vbroadcastss -116(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 193
            vbroadcastss -112(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 193
            vmovaps 18400(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 193
            vmovaps 18432(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 193
            vmovaps 18464(pA), rA2
            #endif

         #endif
         #if KB > 193
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -108(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 194
            vbroadcastss -104(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 194
            vbroadcastss -100(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 194
            vmovaps 18496(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 194
            vmovaps 18528(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 194
            vmovaps 18560(pA), rA2
            #endif

         #endif
         #if KB > 194
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -96(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 195
            vbroadcastss -92(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 195
            vbroadcastss -88(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 195
            vmovaps 18592(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 195
            vmovaps 18624(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 195
            vmovaps 18656(pA), rA2
            #endif

         #endif
         #if KB > 195
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -84(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 196
            vbroadcastss -80(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 196
            vbroadcastss -76(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 196
            vmovaps 18688(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 196
            vmovaps 18720(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 196
            vmovaps 18752(pA), rA2
            #endif

         #endif
         #if KB > 196
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -72(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 197
            vbroadcastss -68(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 197
            vbroadcastss -64(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 197
            vmovaps 18784(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 197
            vmovaps 18816(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 197
            vmovaps 18848(pA), rA2
            #endif

         #endif
         #if KB > 197
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -60(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 198
            vbroadcastss -56(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 198
            vbroadcastss -52(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 198
            vmovaps 18880(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 198
            vmovaps 18912(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 198
            vmovaps 18944(pA), rA2
            #endif

         #endif
         #if KB > 198
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -48(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 199
            vbroadcastss -44(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 199
            vbroadcastss -40(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 199
            vmovaps 18976(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 199
            vmovaps 19008(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 199
            vmovaps 19040(pA), rA2
            #endif

         #endif
         #if KB > 199
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -36(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 200
            vbroadcastss -32(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 200
            vbroadcastss -28(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 200
            vmovaps 19072(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 200
            vmovaps 19104(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 200
            vmovaps 19136(pA), rA2
            #endif

         #endif
         #if KB > 200
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -24(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 201
            vbroadcastss -20(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 201
            vbroadcastss -16(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 201
            vmovaps 19168(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 201
            vmovaps 19200(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 201
            vmovaps 19232(pA), rA2
            #endif

         #endif
         #if KB > 201
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -12(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 202
            vbroadcastss -8(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 202
            vbroadcastss -4(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 202
            vmovaps 19264(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 202
            vmovaps 19296(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 202
            vmovaps 19328(pA), rA2
            #endif

         #endif
         #if KB > 202
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss (pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 203
            vbroadcastss 4(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 203
            vbroadcastss 8(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 203
            vmovaps 19360(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 203
            vmovaps 19392(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 203
            vmovaps 19424(pA), rA2
            #endif

         #endif
         #if KB > 203
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 12(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 204
            vbroadcastss 16(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 204
            vbroadcastss 20(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 204
            vmovaps 19456(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 204
            vmovaps 19488(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 204
            vmovaps 19520(pA), rA2
            #endif

         #endif
         #if KB > 204
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 24(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 205
            vbroadcastss 28(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 205
            vbroadcastss 32(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 205
            vmovaps 19552(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 205
            vmovaps 19584(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 205
            vmovaps 19616(pA), rA2
            #endif

         #endif
         #if KB > 205
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 36(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 206
            vbroadcastss 40(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 206
            vbroadcastss 44(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 206
            vmovaps 19648(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 206
            vmovaps 19680(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 206
            vmovaps 19712(pA), rA2
            #endif

         #endif
         #if KB > 206
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 48(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 207
            vbroadcastss 52(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 207
            vbroadcastss 56(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 207
            vmovaps 19744(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 207
            vmovaps 19776(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 207
            vmovaps 19808(pA), rA2
            #endif

         #endif
         #if KB > 207
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 60(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 208
            vbroadcastss 64(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 208
            vbroadcastss 68(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 208
            vmovaps 19840(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 208
            vmovaps 19872(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 208
            vmovaps 19904(pA), rA2
            #endif

         #endif
         #if KB > 208
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 72(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 209
            vbroadcastss 76(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 209
            vbroadcastss 80(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 209
            vmovaps 19936(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 209
            vmovaps 19968(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 209
            vmovaps 20000(pA), rA2
            #endif

         #endif
         #if KB > 209
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 84(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 210
            vbroadcastss 88(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 210
            vbroadcastss 92(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 210
            vmovaps 20032(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 210
            vmovaps 20064(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 210
            vmovaps 20096(pA), rA2
            #endif

         #endif
         #if KB > 210
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 96(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 211
            vbroadcastss 100(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 211
            vbroadcastss 104(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 211
            vmovaps 20128(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 211
            vmovaps 20160(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 211
            vmovaps 20192(pA), rA2
            #endif

         #endif
         #if KB > 211
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 108(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 212
            vbroadcastss 112(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 212
            vbroadcastss 116(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 212
            vmovaps 20224(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 212
            vmovaps 20256(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 212
            vmovaps 20288(pA), rA2
            #endif

         #endif
         #if KB > 212
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 120(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 213
            vbroadcastss 124(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 213
            vbroadcastss -128(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 213
            vmovaps 20320(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 213
            vmovaps 20352(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 213
            vmovaps 20384(pA), rA2
            #endif

         #endif
         #if KB > 213
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -124(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 214
            vbroadcastss -120(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 214
            vbroadcastss -116(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 214
            vmovaps 20416(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 214
            vmovaps 20448(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 214
            vmovaps 20480(pA), rA2
            #endif

         #endif
         #if KB > 214
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -112(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 215
            vbroadcastss -108(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 215
            vbroadcastss -104(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 215
            vmovaps 20512(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 215
            vmovaps 20544(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 215
            vmovaps 20576(pA), rA2
            #endif

         #endif
         #if KB > 215
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -100(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 216
            vbroadcastss -96(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 216
            vbroadcastss -92(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 216
            vmovaps 20608(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 216
            vmovaps 20640(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 216
            vmovaps 20672(pA), rA2
            #endif

         #endif
         #if KB > 216
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -88(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 217
            vbroadcastss -84(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 217
            vbroadcastss -80(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 217
            vmovaps 20704(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 217
            vmovaps 20736(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 217
            vmovaps 20768(pA), rA2
            #endif

         #endif
         #if KB > 217
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -76(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 218
            vbroadcastss -72(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 218
            vbroadcastss -68(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 218
            vmovaps 20800(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 218
            vmovaps 20832(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 218
            vmovaps 20864(pA), rA2
            #endif

         #endif
         #if KB > 218
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -64(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 219
            vbroadcastss -60(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 219
            vbroadcastss -56(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 219
            vmovaps 20896(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 219
            vmovaps 20928(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 219
            vmovaps 20960(pA), rA2
            #endif

         #endif
         #if KB > 219
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -52(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 220
            vbroadcastss -48(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 220
            vbroadcastss -44(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 220
            vmovaps 20992(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 220
            vmovaps 21024(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 220
            vmovaps 21056(pA), rA2
            #endif

         #endif
         #if KB > 220
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -40(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 221
            vbroadcastss -36(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 221
            vbroadcastss -32(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 221
            vmovaps 21088(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 221
            vmovaps 21120(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 221
            vmovaps 21152(pA), rA2
            #endif

         #endif
         #if KB > 221
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -28(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 222
            vbroadcastss -24(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 222
            vbroadcastss -20(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 222
            vmovaps 21184(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 222
            vmovaps 21216(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 222
            vmovaps 21248(pA), rA2
            #endif

         #endif
         #if KB > 222
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -16(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 223
            vbroadcastss -12(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 223
            vbroadcastss -8(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 223
            vmovaps 21280(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 223
            vmovaps 21312(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 223
            vmovaps 21344(pA), rA2
            #endif

         #endif
         #if KB > 223
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss -4(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 224
            vbroadcastss (pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 224
            vbroadcastss 4(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 224
            vmovaps 21376(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 224
            vmovaps 21408(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 224
            vmovaps 21440(pA), rA2
            #endif

         #endif
         #if KB > 224
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 8(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 225
            vbroadcastss 12(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 225
            vbroadcastss 16(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 225
            vmovaps 21472(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 225
            vmovaps 21504(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 225
            vmovaps 21536(pA), rA2
            #endif

         #endif
         #if KB > 225
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 20(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 226
            vbroadcastss 24(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 226
            vbroadcastss 28(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 226
            vmovaps 21568(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 226
            vmovaps 21600(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 226
            vmovaps 21632(pA), rA2
            #endif

         #endif
         #if KB > 226
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 32(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 227
            vbroadcastss 36(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 227
            vbroadcastss 40(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 227
            vmovaps 21664(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 227
            vmovaps 21696(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 227
            vmovaps 21728(pA), rA2
            #endif

         #endif
         #if KB > 227
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 44(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 228
            vbroadcastss 48(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 228
            vbroadcastss 52(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 228
            vmovaps 21760(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 228
            vmovaps 21792(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 228
            vmovaps 21824(pA), rA2
            #endif

         #endif
         #if KB > 228
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 56(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 229
            vbroadcastss 60(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 229
            vbroadcastss 64(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 229
            vmovaps 21856(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 229
            vmovaps 21888(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 229
            vmovaps 21920(pA), rA2
            #endif

         #endif
         #if KB > 229
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 68(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 230
            vbroadcastss 72(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 230
            vbroadcastss 76(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 230
            vmovaps 21952(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 230
            vmovaps 21984(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 230
            vmovaps 22016(pA), rA2
            #endif

         #endif
         #if KB > 230
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 80(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 231
            vbroadcastss 84(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 231
            vbroadcastss 88(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 231
            vmovaps 22048(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 231
            vmovaps 22080(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 231
            vmovaps 22112(pA), rA2
            #endif

         #endif
         #if KB > 231
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 92(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 232
            vbroadcastss 96(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 232
            vbroadcastss 100(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 232
            vmovaps 22144(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 232
            vmovaps 22176(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 232
            vmovaps 22208(pA), rA2
            #endif

         #endif
         #if KB > 232
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 104(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 233
            vbroadcastss 108(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 233
            vbroadcastss 112(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 233
            vmovaps 22240(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 233
            vmovaps 22272(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 233
            vmovaps 22304(pA), rA2
            #endif

         #endif
         #if KB > 233
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 116(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 234
            vbroadcastss 120(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 234
            vbroadcastss 124(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 234
            vmovaps 22336(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 234
            vmovaps 22368(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 234
            vmovaps 22400(pA), rA2
            #endif

         #endif
         #if KB > 234
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2688(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 235
            vbroadcastss 2692(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 235
            vbroadcastss 2696(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 235
            vmovaps 22432(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 235
            vmovaps 22464(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 235
            vmovaps 22496(pA), rA2
            #endif

         #endif
         #if KB > 235
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2700(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 236
            vbroadcastss 2704(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 236
            vbroadcastss 2708(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 236
            vmovaps 22528(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 236
            vmovaps 22560(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 236
            vmovaps 22592(pA), rA2
            #endif

         #endif
         #if KB > 236
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2712(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 237
            vbroadcastss 2716(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 237
            vbroadcastss 2720(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 237
            vmovaps 22624(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 237
            vmovaps 22656(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 237
            vmovaps 22688(pA), rA2
            #endif

         #endif
         #if KB > 237
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2724(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 238
            vbroadcastss 2728(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 238
            vbroadcastss 2732(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 238
            vmovaps 22720(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 238
            vmovaps 22752(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 238
            vmovaps 22784(pA), rA2
            #endif

         #endif
         #if KB > 238
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2736(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 239
            vbroadcastss 2740(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 239
            vbroadcastss 2744(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 239
            vmovaps 22816(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 239
            vmovaps 22848(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 239
            vmovaps 22880(pA), rA2
            #endif

         #endif
         #if KB > 239
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2748(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 240
            vbroadcastss 2752(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 240
            vbroadcastss 2756(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 240
            vmovaps 22912(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 240
            vmovaps 22944(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 240
            vmovaps 22976(pA), rA2
            #endif

         #endif
         #if KB > 240
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2760(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 241
            vbroadcastss 2764(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 241
            vbroadcastss 2768(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 241
            vmovaps 23008(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 241
            vmovaps 23040(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 241
            vmovaps 23072(pA), rA2
            #endif

         #endif
         #if KB > 241
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2772(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 242
            vbroadcastss 2776(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 242
            vbroadcastss 2780(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 242
            vmovaps 23104(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 242
            vmovaps 23136(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 242
            vmovaps 23168(pA), rA2
            #endif

         #endif
         #if KB > 242
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2784(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 243
            vbroadcastss 2788(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 243
            vbroadcastss 2792(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 243
            vmovaps 23200(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 243
            vmovaps 23232(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 243
            vmovaps 23264(pA), rA2
            #endif

         #endif
         #if KB > 243
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2796(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 244
            vbroadcastss 2800(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 244
            vbroadcastss 2804(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 244
            vmovaps 23296(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 244
            vmovaps 23328(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 244
            vmovaps 23360(pA), rA2
            #endif

         #endif
         #if KB > 244
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2808(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 245
            vbroadcastss 2812(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 245
            vbroadcastss 2816(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 245
            vmovaps 23392(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 245
            vmovaps 23424(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 245
            vmovaps 23456(pA), rA2
            #endif

         #endif
         #if KB > 245
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2820(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 246
            vbroadcastss 2824(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 246
            vbroadcastss 2828(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 246
            vmovaps 23488(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 246
            vmovaps 23520(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 246
            vmovaps 23552(pA), rA2
            #endif

         #endif
         #if KB > 246
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2832(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 247
            vbroadcastss 2836(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 247
            vbroadcastss 2840(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 247
            vmovaps 23584(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 247
            vmovaps 23616(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 247
            vmovaps 23648(pA), rA2
            #endif

         #endif
         #if KB > 247
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2844(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 248
            vbroadcastss 2848(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 248
            vbroadcastss 2852(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 248
            vmovaps 23680(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 248
            vmovaps 23712(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 248
            vmovaps 23744(pA), rA2
            #endif

         #endif
         #if KB > 248
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2856(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 249
            vbroadcastss 2860(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 249
            vbroadcastss 2864(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 249
            vmovaps 23776(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 249
            vmovaps 23808(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 249
            vmovaps 23840(pA), rA2
            #endif

         #endif
         #if KB > 249
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2868(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 250
            vbroadcastss 2872(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 250
            vbroadcastss 2876(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 250
            vmovaps 23872(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 250
            vmovaps 23904(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 250
            vmovaps 23936(pA), rA2
            #endif

         #endif
         #if KB > 250
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2880(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 251
            vbroadcastss 2884(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 251
            vbroadcastss 2888(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 251
            vmovaps 23968(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 251
            vmovaps 24000(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 251
            vmovaps 24032(pA), rA2
            #endif

         #endif
         #if KB > 251
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2892(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 252
            vbroadcastss 2896(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 252
            vbroadcastss 2900(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 252
            vmovaps 24064(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 252
            vmovaps 24096(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 252
            vmovaps 24128(pA), rA2
            #endif

         #endif
         #if KB > 252
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2904(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 253
            vbroadcastss 2908(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 253
            vbroadcastss 2912(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 253
            vmovaps 24160(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 253
            vmovaps 24192(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 253
            vmovaps 24224(pA), rA2
            #endif

         #endif
         #if KB > 253
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2916(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 254
            vbroadcastss 2920(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 254
            vbroadcastss 2924(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 254
            vmovaps 24256(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 254
            vmovaps 24288(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 254
            vmovaps 24320(pA), rA2
            #endif

         #endif
         #if KB > 254
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2928(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 255
            vbroadcastss 2932(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 255
            vbroadcastss 2936(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 255
            vmovaps 24352(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 255
            vmovaps 24384(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 255
            vmovaps 24416(pA), rA2
            #endif

         #endif
         #if KB > 255
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastss 2940(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 256
            vbroadcastss -128(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 256
            vbroadcastss -124(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 256
            vmovaps 24448(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 256
            vmovaps 24480(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 256
            vmovaps 24512(pA), rA2
            #endif

         #endif
         add incPF, pfA
         add incPF, pfB
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
         add $12*3*8, pC        /* pC += MU*NU*sizeof */
         add $KB*3*4, pB        /* pB += K*NU*sizeof */
      sub $1, nnu
      jnz MNLOOP

      mov nnu0, nnu
      mov pB0, pB
      add $KB*12*8, pA          /* pA += KB*MU*size */
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
