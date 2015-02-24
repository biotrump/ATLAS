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
   #define BETCOP vsubpd
#else
   #define BETCOP vaddpd
#endif
#define vmovapd vmovaps
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
            vbroadcastsd -128(pB), rB0
         vmulpd rA0, rB0, rC00
         vmovapd -96(pA), rA1
         vmulpd rA1, rB0, rC10
         vmovapd -64(pA), rA2
         vmulpd rA2, rB0, rC20
            vbroadcastsd -120(pB), rB1
         vmulpd rA0, rB1, rC01
            vbroadcastsd -112(pB), rB2
         vmulpd rA1, rB1, rC11
         #if KB > 1
            vbroadcastsd -104(pB), rB0
         #endif
         vmulpd rA2, rB1, rC21
         #if KB > 1
            vbroadcastsd -96(pB), rB1
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
            vbroadcastsd -88(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefC (pC)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 2
            vbroadcastsd -80(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
               prefC 64(pC)
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 2
            vbroadcastsd -72(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 2
            vmovapd 64(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 2
            vmovapd 96(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 2
            vmovapd -128(pA,i256), rA2
            #endif

         #endif
         #if KB > 2
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -64(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
               prefC 128(pC)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 3
            vbroadcastsd -56(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 3
            vbroadcastsd -48(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 3
            vmovapd -96(pA,i256), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 3
            vmovapd -64(pA,i256), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 3
            vmovapd -32(pA,i256), rA2
            #endif

         #endif
         #if KB > 3
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -40(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefA -128(pfA)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 4
            vbroadcastsd -32(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 4
            vbroadcastsd -24(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 4
            vmovapd (pA,i256), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 4
            vmovapd 32(pA,i256), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 4
            vmovapd 64(pA,i256), rA2
            #endif

         #endif
         #if KB > 4
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -16(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefA -64(pfA)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 5
            vbroadcastsd -8(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 5
            vbroadcastsd (pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 5
            vmovapd 96(pA,i256), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 5
            vmovapd -128(pA,i256,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 5
            vmovapd -96(pA,i256,2), rA2
            #endif

         #endif
         #if KB > 5
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 8(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefA (pfA)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 6
            vbroadcastsd 16(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 6
            vbroadcastsd 24(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 6
            vmovapd -64(pA,i256,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 6
            vmovapd -32(pA,i256,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 6
            vmovapd (pA,i256,2), rA2
            #endif

         #endif
         #if KB > 6
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 32(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefB -128(pfB)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 7
            vbroadcastsd 40(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 7
            vbroadcastsd 48(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 7
            vmovapd 32(pA,i256,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 7
            vmovapd 64(pA,i256,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 7
            vmovapd 96(pA,i256,2), rA2
            #endif

         #endif
         #if KB > 7
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 56(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefB -64(pfB)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 8
            vbroadcastsd 64(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 8
            vbroadcastsd 72(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 8
            vmovapd -128(pA,i768), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 8
            vmovapd -96(pA,i768), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 8
            vmovapd -64(pA,i768), rA2
            #endif

         #endif
         #if KB > 8
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 80(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            prefB (pfB)
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 9
            vbroadcastsd 88(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 9
            vbroadcastsd 96(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 9
            vmovapd -32(pA,i768), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 9
            vmovapd (pA,i768), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 9
            vmovapd 32(pA,i768), rA2
            #endif

         #endif
         #if KB > 9
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 104(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 10
            vbroadcastsd 112(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 10
            vbroadcastsd 120(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 10
            vmovapd 64(pA,i768), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 10
            vmovapd 96(pA,i768), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 10
            vmovapd -128(pA,i256,4), rA2
            #endif

         #endif
         #if KB > 10
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -128(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 11
            vbroadcastsd -120(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 11
            vbroadcastsd -112(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 11
            vmovapd -96(pA,i256,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 11
            vmovapd -64(pA,i256,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 11
            vmovapd -32(pA,i256,4), rA2
            #endif

         #endif
         #if KB > 11
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -104(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 12
            vbroadcastsd -96(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 12
            vbroadcastsd -88(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 12
            vmovapd (pA,i256,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 12
            vmovapd 32(pA,i256,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 12
            vmovapd 64(pA,i256,4), rA2
            #endif

         #endif
         #if KB > 12
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -80(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 13
            vbroadcastsd -72(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 13
            vbroadcastsd -64(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 13
            vmovapd 96(pA,i256,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 13
            vmovapd -128(pA,i1280), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 13
            vmovapd -96(pA,i1280), rA2
            #endif

         #endif
         #if KB > 13
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -56(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 14
            vbroadcastsd -48(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 14
            vbroadcastsd -40(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 14
            vmovapd -64(pA,i1280), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 14
            vmovapd -32(pA,i1280), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 14
            vmovapd (pA,i1280), rA2
            #endif

         #endif
         #if KB > 14
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -32(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 15
            vbroadcastsd -24(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 15
            vbroadcastsd -16(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 15
            vmovapd 32(pA,i1280), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 15
            vmovapd 64(pA,i1280), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 15
            vmovapd 96(pA,i1280), rA2
            #endif

         #endif
         #if KB > 15
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -8(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 16
            vbroadcastsd (pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 16
            vbroadcastsd 8(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 16
            vmovapd -128(pA,i768,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 16
            vmovapd -96(pA,i768,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 16
            vmovapd -64(pA,i768,2), rA2
            #endif

         #endif
         #if KB > 16
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 16(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 17
            vbroadcastsd 24(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 17
            vbroadcastsd 32(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 17
            vmovapd -32(pA,i768,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 17
            vmovapd (pA,i768,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 17
            vmovapd 32(pA,i768,2), rA2
            #endif

         #endif
         #if KB > 17
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 40(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 18
            vbroadcastsd 48(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 18
            vbroadcastsd 56(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 18
            vmovapd 64(pA,i768,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 18
            vmovapd 96(pA,i768,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 18
            vmovapd -128(pA,i1792), rA2
            #endif

         #endif
         #if KB > 18
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 64(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 19
            vbroadcastsd 72(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 19
            vbroadcastsd 80(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 19
            vmovapd -96(pA,i1792), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 19
            vmovapd -64(pA,i1792), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 19
            vmovapd -32(pA,i1792), rA2
            #endif

         #endif
         #if KB > 19
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 88(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 20
            vbroadcastsd 96(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 20
            vbroadcastsd 104(pB,i256), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 20
            vmovapd (pA,i1792), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 20
            vmovapd 32(pA,i1792), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 20
            vmovapd 64(pA,i1792), rA2
            #endif

         #endif
         #if KB > 20
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 112(pB,i256), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 21
            vbroadcastsd 120(pB,i256), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 21
            vbroadcastsd -128(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 21
            vmovapd 96(pA,i1792), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 21
            vmovapd -128(pA,i256,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 21
            vmovapd -96(pA,i256,8), rA2
            #endif

         #endif
         #if KB > 21
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -120(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 22
            vbroadcastsd -112(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 22
            vbroadcastsd -104(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 22
            vmovapd -64(pA,i256,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 22
            vmovapd -32(pA,i256,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 22
            vmovapd (pA,i256,8), rA2
            #endif

         #endif
         #if KB > 22
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -96(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 23
            vbroadcastsd -88(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 23
            vbroadcastsd -80(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 23
            vmovapd 32(pA,i256,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 23
            vmovapd 64(pA,i256,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 23
            vmovapd 96(pA,i256,8), rA2
            #endif

         #endif
         #if KB > 23
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -72(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 24
            vbroadcastsd -64(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 24
            vbroadcastsd -56(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 24
            vmovapd -128(pA,i2304), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 24
            vmovapd -96(pA,i2304), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 24
            vmovapd -64(pA,i2304), rA2
            #endif

         #endif
         #if KB > 24
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -48(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 25
            vbroadcastsd -40(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 25
            vbroadcastsd -32(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 25
            vmovapd -32(pA,i2304), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 25
            vmovapd (pA,i2304), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 25
            vmovapd 32(pA,i2304), rA2
            #endif

         #endif
         #if KB > 25
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -24(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 26
            vbroadcastsd -16(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 26
            vbroadcastsd -8(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 26
            vmovapd 64(pA,i2304), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 26
            vmovapd 96(pA,i2304), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 26
            vmovapd -128(pA,i1280,2), rA2
            #endif

         #endif
         #if KB > 26
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd (pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 27
            vbroadcastsd 8(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 27
            vbroadcastsd 16(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 27
            vmovapd -96(pA,i1280,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 27
            vmovapd -64(pA,i1280,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 27
            vmovapd -32(pA,i1280,2), rA2
            #endif

         #endif
         #if KB > 27
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 24(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 28
            vbroadcastsd 32(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 28
            vbroadcastsd 40(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 28
            vmovapd (pA,i1280,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 28
            vmovapd 32(pA,i1280,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 28
            vmovapd 64(pA,i1280,2), rA2
            #endif

         #endif
         #if KB > 28
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 48(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 29
            vbroadcastsd 56(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 29
            vbroadcastsd 64(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 29
            vmovapd 96(pA,i1280,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 29
            vmovapd 2688(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 29
            vmovapd 2720(pA), rA2
            #endif

         #endif
         #if KB > 29
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 72(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 30
            vbroadcastsd 80(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 30
            vbroadcastsd 88(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 30
            vmovapd 2752(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 30
            vmovapd 2784(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 30
            vmovapd 2816(pA), rA2
            #endif

         #endif
         #if KB > 30
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 96(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 31
            vbroadcastsd 104(pB,i256,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 31
            vbroadcastsd 112(pB,i256,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 31
            vmovapd 2848(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 31
            vmovapd 2880(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 31
            vmovapd 2912(pA), rA2
            #endif

         #endif
         #if KB > 31
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 120(pB,i256,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 32
            vbroadcastsd -128(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 32
            vbroadcastsd -120(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 32
            vmovapd -128(pA,i768,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 32
            vmovapd -96(pA,i768,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 32
            vmovapd -64(pA,i768,4), rA2
            #endif

         #endif
         #if KB > 32
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -112(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 33
            vbroadcastsd -104(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 33
            vbroadcastsd -96(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 33
            vmovapd -32(pA,i768,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 33
            vmovapd (pA,i768,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 33
            vmovapd 32(pA,i768,4), rA2
            #endif

         #endif
         #if KB > 33
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -88(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 34
            vbroadcastsd -80(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 34
            vbroadcastsd -72(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 34
            vmovapd 64(pA,i768,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 34
            vmovapd 96(pA,i768,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 34
            vmovapd 3200(pA), rA2
            #endif

         #endif
         #if KB > 34
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -64(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 35
            vbroadcastsd -56(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 35
            vbroadcastsd -48(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 35
            vmovapd 3232(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 35
            vmovapd 3264(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 35
            vmovapd 3296(pA), rA2
            #endif

         #endif
         #if KB > 35
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -40(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 36
            vbroadcastsd -32(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 36
            vbroadcastsd -24(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 36
            vmovapd 3328(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 36
            vmovapd 3360(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 36
            vmovapd 3392(pA), rA2
            #endif

         #endif
         #if KB > 36
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -16(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 37
            vbroadcastsd -8(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 37
            vbroadcastsd (pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 37
            vmovapd 3424(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 37
            vmovapd -128(pA,i1792,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 37
            vmovapd -96(pA,i1792,2), rA2
            #endif

         #endif
         #if KB > 37
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 8(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 38
            vbroadcastsd 16(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 38
            vbroadcastsd 24(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 38
            vmovapd -64(pA,i1792,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 38
            vmovapd -32(pA,i1792,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 38
            vmovapd (pA,i1792,2), rA2
            #endif

         #endif
         #if KB > 38
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 32(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 39
            vbroadcastsd 40(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 39
            vbroadcastsd 48(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 39
            vmovapd 32(pA,i1792,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 39
            vmovapd 64(pA,i1792,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 39
            vmovapd 96(pA,i1792,2), rA2
            #endif

         #endif
         #if KB > 39
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 56(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 40
            vbroadcastsd 64(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 40
            vbroadcastsd 72(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 40
            vmovapd 3712(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 40
            vmovapd 3744(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 40
            vmovapd 3776(pA), rA2
            #endif

         #endif
         #if KB > 40
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 80(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 41
            vbroadcastsd 88(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 41
            vbroadcastsd 96(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 41
            vmovapd 3808(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 41
            vmovapd 3840(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 41
            vmovapd 3872(pA), rA2
            #endif

         #endif
         #if KB > 41
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 104(pB,i768), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 42
            vbroadcastsd 112(pB,i768), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 42
            vbroadcastsd 120(pB,i768), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 42
            vmovapd 3904(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 42
            vmovapd 3936(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 42
            vmovapd 3968(pA), rA2
            #endif

         #endif
         #if KB > 42
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -128(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 43
            vbroadcastsd -120(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 43
            vbroadcastsd -112(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 43
            vmovapd 4000(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 43
            vmovapd 4032(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 43
            vmovapd 4064(pA), rA2
            #endif

         #endif
         #if KB > 43
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -104(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 44
            vbroadcastsd -96(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 44
            vbroadcastsd -88(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 44
            vmovapd 4096(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 44
            vmovapd 4128(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 44
            vmovapd 4160(pA), rA2
            #endif

         #endif
         #if KB > 44
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -80(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 45
            vbroadcastsd -72(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 45
            vbroadcastsd -64(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 45
            vmovapd 4192(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 45
            vmovapd 4224(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 45
            vmovapd 4256(pA), rA2
            #endif

         #endif
         #if KB > 45
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -56(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 46
            vbroadcastsd -48(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 46
            vbroadcastsd -40(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 46
            vmovapd 4288(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 46
            vmovapd 4320(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 46
            vmovapd 4352(pA), rA2
            #endif

         #endif
         #if KB > 46
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -32(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 47
            vbroadcastsd -24(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 47
            vbroadcastsd -16(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 47
            vmovapd 4384(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 47
            vmovapd 4416(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 47
            vmovapd 4448(pA), rA2
            #endif

         #endif
         #if KB > 47
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -8(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 48
            vbroadcastsd (pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 48
            vbroadcastsd 8(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 48
            vmovapd -128(pA,i2304,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 48
            vmovapd -96(pA,i2304,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 48
            vmovapd -64(pA,i2304,2), rA2
            #endif

         #endif
         #if KB > 48
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 16(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 49
            vbroadcastsd 24(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 49
            vbroadcastsd 32(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 49
            vmovapd -32(pA,i2304,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 49
            vmovapd (pA,i2304,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 49
            vmovapd 32(pA,i2304,2), rA2
            #endif

         #endif
         #if KB > 49
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 40(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 50
            vbroadcastsd 48(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 50
            vbroadcastsd 56(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 50
            vmovapd 64(pA,i2304,2), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 50
            vmovapd 96(pA,i2304,2), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 50
            vmovapd 4736(pA), rA2
            #endif

         #endif
         #if KB > 50
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 64(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 51
            vbroadcastsd 72(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 51
            vbroadcastsd 80(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 51
            vmovapd 4768(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 51
            vmovapd 4800(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 51
            vmovapd 4832(pA), rA2
            #endif

         #endif
         #if KB > 51
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 88(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 52
            vbroadcastsd 96(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 52
            vbroadcastsd 104(pB,i256,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 52
            vmovapd 4864(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 52
            vmovapd 4896(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 52
            vmovapd 4928(pA), rA2
            #endif

         #endif
         #if KB > 52
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 112(pB,i256,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 53
            vbroadcastsd 120(pB,i256,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 53
            vbroadcastsd -128(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 53
            vmovapd 4960(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 53
            vmovapd -128(pA,i1280,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 53
            vmovapd -96(pA,i1280,4), rA2
            #endif

         #endif
         #if KB > 53
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -120(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 54
            vbroadcastsd -112(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 54
            vbroadcastsd -104(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 54
            vmovapd -64(pA,i1280,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 54
            vmovapd -32(pA,i1280,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 54
            vmovapd (pA,i1280,4), rA2
            #endif

         #endif
         #if KB > 54
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -96(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 55
            vbroadcastsd -88(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 55
            vbroadcastsd -80(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 55
            vmovapd 32(pA,i1280,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 55
            vmovapd 64(pA,i1280,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 55
            vmovapd 96(pA,i1280,4), rA2
            #endif

         #endif
         #if KB > 55
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -72(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 56
            vbroadcastsd -64(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 56
            vbroadcastsd -56(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 56
            vmovapd 5248(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 56
            vmovapd 5280(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 56
            vmovapd 5312(pA), rA2
            #endif

         #endif
         #if KB > 56
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -48(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 57
            vbroadcastsd -40(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 57
            vbroadcastsd -32(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 57
            vmovapd 5344(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 57
            vmovapd 5376(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 57
            vmovapd 5408(pA), rA2
            #endif

         #endif
         #if KB > 57
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -24(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 58
            vbroadcastsd -16(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 58
            vbroadcastsd -8(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 58
            vmovapd 5440(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 58
            vmovapd 5472(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 58
            vmovapd 5504(pA), rA2
            #endif

         #endif
         #if KB > 58
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd (pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 59
            vbroadcastsd 8(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 59
            vbroadcastsd 16(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 59
            vmovapd 5536(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 59
            vmovapd 5568(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 59
            vmovapd 5600(pA), rA2
            #endif

         #endif
         #if KB > 59
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 24(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 60
            vbroadcastsd 32(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 60
            vbroadcastsd 40(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 60
            vmovapd 5632(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 60
            vmovapd 5664(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 60
            vmovapd 5696(pA), rA2
            #endif

         #endif
         #if KB > 60
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 48(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 61
            vbroadcastsd 56(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 61
            vbroadcastsd 64(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 61
            vmovapd 5728(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 61
            vmovapd 5760(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 61
            vmovapd 5792(pA), rA2
            #endif

         #endif
         #if KB > 61
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 72(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 62
            vbroadcastsd 80(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 62
            vbroadcastsd 88(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 62
            vmovapd 5824(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 62
            vmovapd 5856(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 62
            vmovapd 5888(pA), rA2
            #endif

         #endif
         #if KB > 62
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 96(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 63
            vbroadcastsd 104(pB,i1280), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 63
            vbroadcastsd 112(pB,i1280), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 63
            vmovapd 5920(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 63
            vmovapd 5952(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 63
            vmovapd 5984(pA), rA2
            #endif

         #endif
         #if KB > 63
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 120(pB,i1280), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 64
            vbroadcastsd -128(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 64
            vbroadcastsd -120(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 64
            vmovapd -128(pA,i768,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 64
            vmovapd -96(pA,i768,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 64
            vmovapd -64(pA,i768,8), rA2
            #endif

         #endif
         #if KB > 64
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -112(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 65
            vbroadcastsd -104(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 65
            vbroadcastsd -96(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 65
            vmovapd -32(pA,i768,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 65
            vmovapd (pA,i768,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 65
            vmovapd 32(pA,i768,8), rA2
            #endif

         #endif
         #if KB > 65
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -88(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 66
            vbroadcastsd -80(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 66
            vbroadcastsd -72(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 66
            vmovapd 64(pA,i768,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 66
            vmovapd 96(pA,i768,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 66
            vmovapd 6272(pA), rA2
            #endif

         #endif
         #if KB > 66
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -64(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 67
            vbroadcastsd -56(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 67
            vbroadcastsd -48(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 67
            vmovapd 6304(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 67
            vmovapd 6336(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 67
            vmovapd 6368(pA), rA2
            #endif

         #endif
         #if KB > 67
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -40(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 68
            vbroadcastsd -32(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 68
            vbroadcastsd -24(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 68
            vmovapd 6400(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 68
            vmovapd 6432(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 68
            vmovapd 6464(pA), rA2
            #endif

         #endif
         #if KB > 68
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -16(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 69
            vbroadcastsd -8(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 69
            vbroadcastsd (pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 69
            vmovapd 6496(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 69
            vmovapd 6528(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 69
            vmovapd 6560(pA), rA2
            #endif

         #endif
         #if KB > 69
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 8(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 70
            vbroadcastsd 16(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 70
            vbroadcastsd 24(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 70
            vmovapd 6592(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 70
            vmovapd 6624(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 70
            vmovapd 6656(pA), rA2
            #endif

         #endif
         #if KB > 70
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 32(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 71
            vbroadcastsd 40(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 71
            vbroadcastsd 48(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 71
            vmovapd 6688(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 71
            vmovapd 6720(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 71
            vmovapd 6752(pA), rA2
            #endif

         #endif
         #if KB > 71
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 56(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 72
            vbroadcastsd 64(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 72
            vbroadcastsd 72(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 72
            vmovapd 6784(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 72
            vmovapd 6816(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 72
            vmovapd 6848(pA), rA2
            #endif

         #endif
         #if KB > 72
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 80(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 73
            vbroadcastsd 88(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 73
            vbroadcastsd 96(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 73
            vmovapd 6880(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 73
            vmovapd 6912(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 73
            vmovapd 6944(pA), rA2
            #endif

         #endif
         #if KB > 73
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 104(pB,i768,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 74
            vbroadcastsd 112(pB,i768,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 74
            vbroadcastsd 120(pB,i768,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 74
            vmovapd 6976(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 74
            vmovapd 7008(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 74
            vmovapd -128(pA,i1792,4), rA2
            #endif

         #endif
         #if KB > 74
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -128(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 75
            vbroadcastsd -120(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 75
            vbroadcastsd -112(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 75
            vmovapd -96(pA,i1792,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 75
            vmovapd -64(pA,i1792,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 75
            vmovapd -32(pA,i1792,4), rA2
            #endif

         #endif
         #if KB > 75
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -104(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 76
            vbroadcastsd -96(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 76
            vbroadcastsd -88(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 76
            vmovapd (pA,i1792,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 76
            vmovapd 32(pA,i1792,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 76
            vmovapd 64(pA,i1792,4), rA2
            #endif

         #endif
         #if KB > 76
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -80(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 77
            vbroadcastsd -72(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 77
            vbroadcastsd -64(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 77
            vmovapd 96(pA,i1792,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 77
            vmovapd 7296(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 77
            vmovapd 7328(pA), rA2
            #endif

         #endif
         #if KB > 77
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -56(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 78
            vbroadcastsd -48(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 78
            vbroadcastsd -40(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 78
            vmovapd 7360(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 78
            vmovapd 7392(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 78
            vmovapd 7424(pA), rA2
            #endif

         #endif
         #if KB > 78
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -32(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 79
            vbroadcastsd -24(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 79
            vbroadcastsd -16(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 79
            vmovapd 7456(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 79
            vmovapd 7488(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 79
            vmovapd 7520(pA), rA2
            #endif

         #endif
         #if KB > 79
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -8(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 80
            vbroadcastsd (pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 80
            vbroadcastsd 8(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 80
            vmovapd 7552(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 80
            vmovapd 7584(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 80
            vmovapd 7616(pA), rA2
            #endif

         #endif
         #if KB > 80
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 16(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 81
            vbroadcastsd 24(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 81
            vbroadcastsd 32(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 81
            vmovapd 7648(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 81
            vmovapd 7680(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 81
            vmovapd 7712(pA), rA2
            #endif

         #endif
         #if KB > 81
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 40(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 82
            vbroadcastsd 48(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 82
            vbroadcastsd 56(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 82
            vmovapd 7744(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 82
            vmovapd 7776(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 82
            vmovapd 7808(pA), rA2
            #endif

         #endif
         #if KB > 82
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 64(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 83
            vbroadcastsd 72(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 83
            vbroadcastsd 80(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 83
            vmovapd 7840(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 83
            vmovapd 7872(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 83
            vmovapd 7904(pA), rA2
            #endif

         #endif
         #if KB > 83
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 88(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 84
            vbroadcastsd 96(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 84
            vbroadcastsd 104(pB,i1792), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 84
            vmovapd 7936(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 84
            vmovapd 7968(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 84
            vmovapd 8000(pA), rA2
            #endif

         #endif
         #if KB > 84
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 112(pB,i1792), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 85
            vbroadcastsd 120(pB,i1792), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 85
            vbroadcastsd -128(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 85
            vmovapd 8032(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 85
            vmovapd 8064(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 85
            vmovapd 8096(pA), rA2
            #endif

         #endif
         #if KB > 85
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -120(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 86
            vbroadcastsd -112(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 86
            vbroadcastsd -104(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 86
            vmovapd 8128(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 86
            vmovapd 8160(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 86
            vmovapd 8192(pA), rA2
            #endif

         #endif
         #if KB > 86
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -96(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 87
            vbroadcastsd -88(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 87
            vbroadcastsd -80(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 87
            vmovapd 8224(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 87
            vmovapd 8256(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 87
            vmovapd 8288(pA), rA2
            #endif

         #endif
         #if KB > 87
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -72(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 88
            vbroadcastsd -64(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 88
            vbroadcastsd -56(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 88
            vmovapd 8320(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 88
            vmovapd 8352(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 88
            vmovapd 8384(pA), rA2
            #endif

         #endif
         #if KB > 88
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -48(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 89
            vbroadcastsd -40(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 89
            vbroadcastsd -32(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 89
            vmovapd 8416(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 89
            vmovapd 8448(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 89
            vmovapd 8480(pA), rA2
            #endif

         #endif
         #if KB > 89
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -24(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 90
            vbroadcastsd -16(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 90
            vbroadcastsd -8(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 90
            vmovapd 8512(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 90
            vmovapd 8544(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 90
            vmovapd 8576(pA), rA2
            #endif

         #endif
         #if KB > 90
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd (pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 91
            vbroadcastsd 8(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 91
            vbroadcastsd 16(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 91
            vmovapd 8608(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 91
            vmovapd 8640(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 91
            vmovapd 8672(pA), rA2
            #endif

         #endif
         #if KB > 91
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 24(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 92
            vbroadcastsd 32(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 92
            vbroadcastsd 40(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 92
            vmovapd 8704(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 92
            vmovapd 8736(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 92
            vmovapd 8768(pA), rA2
            #endif

         #endif
         #if KB > 92
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 48(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 93
            vbroadcastsd 56(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 93
            vbroadcastsd 64(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 93
            vmovapd 8800(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 93
            vmovapd 8832(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 93
            vmovapd 8864(pA), rA2
            #endif

         #endif
         #if KB > 93
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 72(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 94
            vbroadcastsd 80(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 94
            vbroadcastsd 88(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 94
            vmovapd 8896(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 94
            vmovapd 8928(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 94
            vmovapd 8960(pA), rA2
            #endif

         #endif
         #if KB > 94
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 96(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 95
            vbroadcastsd 104(pB,i256,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 95
            vbroadcastsd 112(pB,i256,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 95
            vmovapd 8992(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 95
            vmovapd 9024(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 95
            vmovapd 9056(pA), rA2
            #endif

         #endif
         #if KB > 95
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 120(pB,i256,8), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 96
            vbroadcastsd -128(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 96
            vbroadcastsd -120(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 96
            vmovapd -128(pA,i2304,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 96
            vmovapd -96(pA,i2304,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 96
            vmovapd -64(pA,i2304,4), rA2
            #endif

         #endif
         #if KB > 96
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -112(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 97
            vbroadcastsd -104(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 97
            vbroadcastsd -96(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 97
            vmovapd -32(pA,i2304,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 97
            vmovapd (pA,i2304,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 97
            vmovapd 32(pA,i2304,4), rA2
            #endif

         #endif
         #if KB > 97
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -88(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 98
            vbroadcastsd -80(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 98
            vbroadcastsd -72(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 98
            vmovapd 64(pA,i2304,4), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 98
            vmovapd 96(pA,i2304,4), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 98
            vmovapd 9344(pA), rA2
            #endif

         #endif
         #if KB > 98
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -64(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 99
            vbroadcastsd -56(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 99
            vbroadcastsd -48(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 99
            vmovapd 9376(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 99
            vmovapd 9408(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 99
            vmovapd 9440(pA), rA2
            #endif

         #endif
         #if KB > 99
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -40(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 100
            vbroadcastsd -32(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 100
            vbroadcastsd -24(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 100
            vmovapd 9472(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 100
            vmovapd 9504(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 100
            vmovapd 9536(pA), rA2
            #endif

         #endif
         #if KB > 100
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -16(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 101
            vbroadcastsd -8(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 101
            vbroadcastsd (pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 101
            vmovapd 9568(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 101
            vmovapd 9600(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 101
            vmovapd 9632(pA), rA2
            #endif

         #endif
         #if KB > 101
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 8(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 102
            vbroadcastsd 16(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 102
            vbroadcastsd 24(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 102
            vmovapd 9664(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 102
            vmovapd 9696(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 102
            vmovapd 9728(pA), rA2
            #endif

         #endif
         #if KB > 102
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 32(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 103
            vbroadcastsd 40(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 103
            vbroadcastsd 48(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 103
            vmovapd 9760(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 103
            vmovapd 9792(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 103
            vmovapd 9824(pA), rA2
            #endif

         #endif
         #if KB > 103
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 56(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 104
            vbroadcastsd 64(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 104
            vbroadcastsd 72(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 104
            vmovapd 9856(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 104
            vmovapd 9888(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 104
            vmovapd 9920(pA), rA2
            #endif

         #endif
         #if KB > 104
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 80(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 105
            vbroadcastsd 88(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 105
            vbroadcastsd 96(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 105
            vmovapd 9952(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 105
            vmovapd 9984(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 105
            vmovapd 10016(pA), rA2
            #endif

         #endif
         #if KB > 105
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 104(pB,i2304), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 106
            vbroadcastsd 112(pB,i2304), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 106
            vbroadcastsd 120(pB,i2304), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 106
            vmovapd 10048(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 106
            vmovapd 10080(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 106
            vmovapd -128(pA,i1280,8), rA2
            #endif

         #endif
         #if KB > 106
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -128(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 107
            vbroadcastsd -120(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 107
            vbroadcastsd -112(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 107
            vmovapd -96(pA,i1280,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 107
            vmovapd -64(pA,i1280,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 107
            vmovapd -32(pA,i1280,8), rA2
            #endif

         #endif
         #if KB > 107
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -104(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 108
            vbroadcastsd -96(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 108
            vbroadcastsd -88(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 108
            vmovapd (pA,i1280,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 108
            vmovapd 32(pA,i1280,8), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 108
            vmovapd 64(pA,i1280,8), rA2
            #endif

         #endif
         #if KB > 108
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -80(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 109
            vbroadcastsd -72(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 109
            vbroadcastsd -64(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 109
            vmovapd 96(pA,i1280,8), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 109
            vmovapd 10368(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 109
            vmovapd 10400(pA), rA2
            #endif

         #endif
         #if KB > 109
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -56(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 110
            vbroadcastsd -48(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 110
            vbroadcastsd -40(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 110
            vmovapd 10432(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 110
            vmovapd 10464(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 110
            vmovapd 10496(pA), rA2
            #endif

         #endif
         #if KB > 110
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -32(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 111
            vbroadcastsd -24(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 111
            vbroadcastsd -16(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 111
            vmovapd 10528(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 111
            vmovapd 10560(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 111
            vmovapd 10592(pA), rA2
            #endif

         #endif
         #if KB > 111
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -8(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 112
            vbroadcastsd (pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 112
            vbroadcastsd 8(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 112
            vmovapd 10624(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 112
            vmovapd 10656(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 112
            vmovapd 10688(pA), rA2
            #endif

         #endif
         #if KB > 112
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 16(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 113
            vbroadcastsd 24(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 113
            vbroadcastsd 32(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 113
            vmovapd 10720(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 113
            vmovapd 10752(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 113
            vmovapd 10784(pA), rA2
            #endif

         #endif
         #if KB > 113
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 40(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 114
            vbroadcastsd 48(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 114
            vbroadcastsd 56(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 114
            vmovapd 10816(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 114
            vmovapd 10848(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 114
            vmovapd 10880(pA), rA2
            #endif

         #endif
         #if KB > 114
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 64(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 115
            vbroadcastsd 72(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 115
            vbroadcastsd 80(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 115
            vmovapd 10912(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 115
            vmovapd 10944(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 115
            vmovapd 10976(pA), rA2
            #endif

         #endif
         #if KB > 115
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 88(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 116
            vbroadcastsd 96(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 116
            vbroadcastsd 104(pB,i1280,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 116
            vmovapd 11008(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 116
            vmovapd 11040(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 116
            vmovapd 11072(pA), rA2
            #endif

         #endif
         #if KB > 116
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 112(pB,i1280,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 117
            vbroadcastsd 120(pB,i1280,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 117
            vbroadcastsd 2688(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 117
            vmovapd 11104(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 117
            vmovapd 11136(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 117
            vmovapd 11168(pA), rA2
            #endif

         #endif
         #if KB > 117
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2696(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 118
            vbroadcastsd 2704(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 118
            vbroadcastsd 2712(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 118
            vmovapd 11200(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 118
            vmovapd 11232(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 118
            vmovapd 11264(pA), rA2
            #endif

         #endif
         #if KB > 118
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2720(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 119
            vbroadcastsd 2728(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 119
            vbroadcastsd 2736(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 119
            vmovapd 11296(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 119
            vmovapd 11328(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 119
            vmovapd 11360(pA), rA2
            #endif

         #endif
         #if KB > 119
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2744(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 120
            vbroadcastsd 2752(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 120
            vbroadcastsd 2760(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 120
            vmovapd 11392(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 120
            vmovapd 11424(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 120
            vmovapd 11456(pA), rA2
            #endif

         #endif
         #if KB > 120
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2768(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 121
            vbroadcastsd 2776(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 121
            vbroadcastsd 2784(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 121
            vmovapd 11488(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 121
            vmovapd 11520(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 121
            vmovapd 11552(pA), rA2
            #endif

         #endif
         #if KB > 121
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2792(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 122
            vbroadcastsd 2800(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 122
            vbroadcastsd 2808(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 122
            vmovapd 11584(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 122
            vmovapd 11616(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 122
            vmovapd 11648(pA), rA2
            #endif

         #endif
         #if KB > 122
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2816(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 123
            vbroadcastsd 2824(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 123
            vbroadcastsd 2832(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 123
            vmovapd 11680(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 123
            vmovapd 11712(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 123
            vmovapd 11744(pA), rA2
            #endif

         #endif
         #if KB > 123
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2840(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 124
            vbroadcastsd 2848(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 124
            vbroadcastsd 2856(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 124
            vmovapd 11776(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 124
            vmovapd 11808(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 124
            vmovapd 11840(pA), rA2
            #endif

         #endif
         #if KB > 124
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2864(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 125
            vbroadcastsd 2872(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 125
            vbroadcastsd 2880(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 125
            vmovapd 11872(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 125
            vmovapd 11904(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 125
            vmovapd 11936(pA), rA2
            #endif

         #endif
         #if KB > 125
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2888(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 126
            vbroadcastsd 2896(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 126
            vbroadcastsd 2904(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 126
            vmovapd 11968(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 126
            vmovapd 12000(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 126
            vmovapd 12032(pA), rA2
            #endif

         #endif
         #if KB > 126
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2912(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 127
            vbroadcastsd 2920(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 127
            vbroadcastsd 2928(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 127
            vmovapd 12064(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 127
            vmovapd 12096(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 127
            vmovapd 12128(pA), rA2
            #endif

         #endif
         #if KB > 127
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 2936(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 128
            vbroadcastsd -128(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 128
            vbroadcastsd -120(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 128
            vmovapd 12160(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 128
            vmovapd 12192(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 128
            vmovapd 12224(pA), rA2
            #endif

         #endif
         #if KB > 128
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -112(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 129
            vbroadcastsd -104(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 129
            vbroadcastsd -96(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 129
            vmovapd 12256(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 129
            vmovapd 12288(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 129
            vmovapd 12320(pA), rA2
            #endif

         #endif
         #if KB > 129
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -88(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 130
            vbroadcastsd -80(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 130
            vbroadcastsd -72(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 130
            vmovapd 12352(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 130
            vmovapd 12384(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 130
            vmovapd 12416(pA), rA2
            #endif

         #endif
         #if KB > 130
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -64(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 131
            vbroadcastsd -56(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 131
            vbroadcastsd -48(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 131
            vmovapd 12448(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 131
            vmovapd 12480(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 131
            vmovapd 12512(pA), rA2
            #endif

         #endif
         #if KB > 131
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -40(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 132
            vbroadcastsd -32(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 132
            vbroadcastsd -24(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 132
            vmovapd 12544(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 132
            vmovapd 12576(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 132
            vmovapd 12608(pA), rA2
            #endif

         #endif
         #if KB > 132
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -16(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 133
            vbroadcastsd -8(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 133
            vbroadcastsd (pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 133
            vmovapd 12640(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 133
            vmovapd 12672(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 133
            vmovapd 12704(pA), rA2
            #endif

         #endif
         #if KB > 133
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 8(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 134
            vbroadcastsd 16(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 134
            vbroadcastsd 24(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 134
            vmovapd 12736(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 134
            vmovapd 12768(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 134
            vmovapd 12800(pA), rA2
            #endif

         #endif
         #if KB > 134
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 32(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 135
            vbroadcastsd 40(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 135
            vbroadcastsd 48(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 135
            vmovapd 12832(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 135
            vmovapd 12864(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 135
            vmovapd 12896(pA), rA2
            #endif

         #endif
         #if KB > 135
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 56(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 136
            vbroadcastsd 64(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 136
            vbroadcastsd 72(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 136
            vmovapd 12928(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 136
            vmovapd 12960(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 136
            vmovapd 12992(pA), rA2
            #endif

         #endif
         #if KB > 136
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 80(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 137
            vbroadcastsd 88(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 137
            vbroadcastsd 96(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 137
            vmovapd 13024(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 137
            vmovapd 13056(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 137
            vmovapd 13088(pA), rA2
            #endif

         #endif
         #if KB > 137
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 104(pB,i768,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 138
            vbroadcastsd 112(pB,i768,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 138
            vbroadcastsd 120(pB,i768,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 138
            vmovapd 13120(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 138
            vmovapd 13152(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 138
            vmovapd 13184(pA), rA2
            #endif

         #endif
         #if KB > 138
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3200(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 139
            vbroadcastsd 3208(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 139
            vbroadcastsd 3216(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 139
            vmovapd 13216(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 139
            vmovapd 13248(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 139
            vmovapd 13280(pA), rA2
            #endif

         #endif
         #if KB > 139
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3224(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 140
            vbroadcastsd 3232(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 140
            vbroadcastsd 3240(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 140
            vmovapd 13312(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 140
            vmovapd 13344(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 140
            vmovapd 13376(pA), rA2
            #endif

         #endif
         #if KB > 140
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3248(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 141
            vbroadcastsd 3256(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 141
            vbroadcastsd 3264(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 141
            vmovapd 13408(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 141
            vmovapd 13440(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 141
            vmovapd 13472(pA), rA2
            #endif

         #endif
         #if KB > 141
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3272(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 142
            vbroadcastsd 3280(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 142
            vbroadcastsd 3288(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 142
            vmovapd 13504(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 142
            vmovapd 13536(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 142
            vmovapd 13568(pA), rA2
            #endif

         #endif
         #if KB > 142
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3296(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 143
            vbroadcastsd 3304(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 143
            vbroadcastsd 3312(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 143
            vmovapd 13600(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 143
            vmovapd 13632(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 143
            vmovapd 13664(pA), rA2
            #endif

         #endif
         #if KB > 143
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3320(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 144
            vbroadcastsd 3328(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 144
            vbroadcastsd 3336(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 144
            vmovapd 13696(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 144
            vmovapd 13728(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 144
            vmovapd 13760(pA), rA2
            #endif

         #endif
         #if KB > 144
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3344(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 145
            vbroadcastsd 3352(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 145
            vbroadcastsd 3360(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 145
            vmovapd 13792(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 145
            vmovapd 13824(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 145
            vmovapd 13856(pA), rA2
            #endif

         #endif
         #if KB > 145
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3368(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 146
            vbroadcastsd 3376(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 146
            vbroadcastsd 3384(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 146
            vmovapd 13888(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 146
            vmovapd 13920(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 146
            vmovapd 13952(pA), rA2
            #endif

         #endif
         #if KB > 146
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3392(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 147
            vbroadcastsd 3400(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 147
            vbroadcastsd 3408(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 147
            vmovapd 13984(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 147
            vmovapd 14016(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 147
            vmovapd 14048(pA), rA2
            #endif

         #endif
         #if KB > 147
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3416(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 148
            vbroadcastsd 3424(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 148
            vbroadcastsd 3432(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 148
            vmovapd 14080(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 148
            vmovapd 14112(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 148
            vmovapd 14144(pA), rA2
            #endif

         #endif
         #if KB > 148
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3440(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 149
            vbroadcastsd 3448(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 149
            vbroadcastsd -128(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 149
            vmovapd 14176(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 149
            vmovapd 14208(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 149
            vmovapd 14240(pA), rA2
            #endif

         #endif
         #if KB > 149
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -120(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 150
            vbroadcastsd -112(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 150
            vbroadcastsd -104(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 150
            vmovapd 14272(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 150
            vmovapd 14304(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 150
            vmovapd 14336(pA), rA2
            #endif

         #endif
         #if KB > 150
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -96(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 151
            vbroadcastsd -88(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 151
            vbroadcastsd -80(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 151
            vmovapd 14368(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 151
            vmovapd 14400(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 151
            vmovapd 14432(pA), rA2
            #endif

         #endif
         #if KB > 151
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -72(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 152
            vbroadcastsd -64(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 152
            vbroadcastsd -56(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 152
            vmovapd 14464(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 152
            vmovapd 14496(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 152
            vmovapd 14528(pA), rA2
            #endif

         #endif
         #if KB > 152
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -48(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 153
            vbroadcastsd -40(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 153
            vbroadcastsd -32(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 153
            vmovapd 14560(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 153
            vmovapd 14592(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 153
            vmovapd 14624(pA), rA2
            #endif

         #endif
         #if KB > 153
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -24(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 154
            vbroadcastsd -16(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 154
            vbroadcastsd -8(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 154
            vmovapd 14656(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 154
            vmovapd 14688(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 154
            vmovapd 14720(pA), rA2
            #endif

         #endif
         #if KB > 154
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd (pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 155
            vbroadcastsd 8(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 155
            vbroadcastsd 16(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 155
            vmovapd 14752(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 155
            vmovapd 14784(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 155
            vmovapd 14816(pA), rA2
            #endif

         #endif
         #if KB > 155
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 24(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 156
            vbroadcastsd 32(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 156
            vbroadcastsd 40(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 156
            vmovapd 14848(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 156
            vmovapd 14880(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 156
            vmovapd 14912(pA), rA2
            #endif

         #endif
         #if KB > 156
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 48(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 157
            vbroadcastsd 56(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 157
            vbroadcastsd 64(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 157
            vmovapd 14944(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 157
            vmovapd 14976(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 157
            vmovapd 15008(pA), rA2
            #endif

         #endif
         #if KB > 157
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 72(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 158
            vbroadcastsd 80(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 158
            vbroadcastsd 88(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 158
            vmovapd 15040(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 158
            vmovapd 15072(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 158
            vmovapd 15104(pA), rA2
            #endif

         #endif
         #if KB > 158
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 96(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 159
            vbroadcastsd 104(pB,i1792,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 159
            vbroadcastsd 112(pB,i1792,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 159
            vmovapd 15136(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 159
            vmovapd 15168(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 159
            vmovapd 15200(pA), rA2
            #endif

         #endif
         #if KB > 159
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 120(pB,i1792,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 160
            vbroadcastsd 3712(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 160
            vbroadcastsd 3720(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 160
            vmovapd 15232(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 160
            vmovapd 15264(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 160
            vmovapd 15296(pA), rA2
            #endif

         #endif
         #if KB > 160
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3728(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 161
            vbroadcastsd 3736(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 161
            vbroadcastsd 3744(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 161
            vmovapd 15328(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 161
            vmovapd 15360(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 161
            vmovapd 15392(pA), rA2
            #endif

         #endif
         #if KB > 161
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3752(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 162
            vbroadcastsd 3760(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 162
            vbroadcastsd 3768(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 162
            vmovapd 15424(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 162
            vmovapd 15456(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 162
            vmovapd 15488(pA), rA2
            #endif

         #endif
         #if KB > 162
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3776(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 163
            vbroadcastsd 3784(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 163
            vbroadcastsd 3792(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 163
            vmovapd 15520(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 163
            vmovapd 15552(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 163
            vmovapd 15584(pA), rA2
            #endif

         #endif
         #if KB > 163
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3800(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 164
            vbroadcastsd 3808(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 164
            vbroadcastsd 3816(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 164
            vmovapd 15616(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 164
            vmovapd 15648(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 164
            vmovapd 15680(pA), rA2
            #endif

         #endif
         #if KB > 164
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3824(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 165
            vbroadcastsd 3832(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 165
            vbroadcastsd 3840(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 165
            vmovapd 15712(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 165
            vmovapd 15744(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 165
            vmovapd 15776(pA), rA2
            #endif

         #endif
         #if KB > 165
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3848(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 166
            vbroadcastsd 3856(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 166
            vbroadcastsd 3864(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 166
            vmovapd 15808(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 166
            vmovapd 15840(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 166
            vmovapd 15872(pA), rA2
            #endif

         #endif
         #if KB > 166
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3872(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 167
            vbroadcastsd 3880(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 167
            vbroadcastsd 3888(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 167
            vmovapd 15904(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 167
            vmovapd 15936(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 167
            vmovapd 15968(pA), rA2
            #endif

         #endif
         #if KB > 167
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3896(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 168
            vbroadcastsd 3904(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 168
            vbroadcastsd 3912(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 168
            vmovapd 16000(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 168
            vmovapd 16032(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 168
            vmovapd 16064(pA), rA2
            #endif

         #endif
         #if KB > 168
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3920(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 169
            vbroadcastsd 3928(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 169
            vbroadcastsd 3936(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 169
            vmovapd 16096(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 169
            vmovapd 16128(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 169
            vmovapd 16160(pA), rA2
            #endif

         #endif
         #if KB > 169
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3944(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 170
            vbroadcastsd 3952(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 170
            vbroadcastsd 3960(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 170
            vmovapd 16192(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 170
            vmovapd 16224(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 170
            vmovapd 16256(pA), rA2
            #endif

         #endif
         #if KB > 170
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3968(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 171
            vbroadcastsd 3976(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 171
            vbroadcastsd 3984(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 171
            vmovapd 16288(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 171
            vmovapd 16320(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 171
            vmovapd 16352(pA), rA2
            #endif

         #endif
         #if KB > 171
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 3992(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 172
            vbroadcastsd 4000(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 172
            vbroadcastsd 4008(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 172
            vmovapd 16384(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 172
            vmovapd 16416(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 172
            vmovapd 16448(pA), rA2
            #endif

         #endif
         #if KB > 172
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4016(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 173
            vbroadcastsd 4024(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 173
            vbroadcastsd 4032(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 173
            vmovapd 16480(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 173
            vmovapd 16512(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 173
            vmovapd 16544(pA), rA2
            #endif

         #endif
         #if KB > 173
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4040(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 174
            vbroadcastsd 4048(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 174
            vbroadcastsd 4056(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 174
            vmovapd 16576(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 174
            vmovapd 16608(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 174
            vmovapd 16640(pA), rA2
            #endif

         #endif
         #if KB > 174
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4064(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 175
            vbroadcastsd 4072(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 175
            vbroadcastsd 4080(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 175
            vmovapd 16672(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 175
            vmovapd 16704(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 175
            vmovapd 16736(pA), rA2
            #endif

         #endif
         #if KB > 175
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4088(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 176
            vbroadcastsd 4096(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 176
            vbroadcastsd 4104(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 176
            vmovapd 16768(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 176
            vmovapd 16800(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 176
            vmovapd 16832(pA), rA2
            #endif

         #endif
         #if KB > 176
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4112(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 177
            vbroadcastsd 4120(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 177
            vbroadcastsd 4128(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 177
            vmovapd 16864(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 177
            vmovapd 16896(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 177
            vmovapd 16928(pA), rA2
            #endif

         #endif
         #if KB > 177
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4136(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 178
            vbroadcastsd 4144(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 178
            vbroadcastsd 4152(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 178
            vmovapd 16960(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 178
            vmovapd 16992(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 178
            vmovapd 17024(pA), rA2
            #endif

         #endif
         #if KB > 178
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4160(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 179
            vbroadcastsd 4168(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 179
            vbroadcastsd 4176(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 179
            vmovapd 17056(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 179
            vmovapd 17088(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 179
            vmovapd 17120(pA), rA2
            #endif

         #endif
         #if KB > 179
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4184(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 180
            vbroadcastsd 4192(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 180
            vbroadcastsd 4200(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 180
            vmovapd 17152(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 180
            vmovapd 17184(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 180
            vmovapd 17216(pA), rA2
            #endif

         #endif
         #if KB > 180
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4208(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 181
            vbroadcastsd 4216(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 181
            vbroadcastsd 4224(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 181
            vmovapd 17248(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 181
            vmovapd 17280(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 181
            vmovapd 17312(pA), rA2
            #endif

         #endif
         #if KB > 181
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4232(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 182
            vbroadcastsd 4240(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 182
            vbroadcastsd 4248(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 182
            vmovapd 17344(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 182
            vmovapd 17376(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 182
            vmovapd 17408(pA), rA2
            #endif

         #endif
         #if KB > 182
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4256(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 183
            vbroadcastsd 4264(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 183
            vbroadcastsd 4272(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 183
            vmovapd 17440(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 183
            vmovapd 17472(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 183
            vmovapd 17504(pA), rA2
            #endif

         #endif
         #if KB > 183
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4280(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 184
            vbroadcastsd 4288(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 184
            vbroadcastsd 4296(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 184
            vmovapd 17536(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 184
            vmovapd 17568(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 184
            vmovapd 17600(pA), rA2
            #endif

         #endif
         #if KB > 184
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4304(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 185
            vbroadcastsd 4312(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 185
            vbroadcastsd 4320(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 185
            vmovapd 17632(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 185
            vmovapd 17664(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 185
            vmovapd 17696(pA), rA2
            #endif

         #endif
         #if KB > 185
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4328(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 186
            vbroadcastsd 4336(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 186
            vbroadcastsd 4344(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 186
            vmovapd 17728(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 186
            vmovapd 17760(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 186
            vmovapd 17792(pA), rA2
            #endif

         #endif
         #if KB > 186
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4352(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 187
            vbroadcastsd 4360(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 187
            vbroadcastsd 4368(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 187
            vmovapd 17824(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 187
            vmovapd 17856(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 187
            vmovapd 17888(pA), rA2
            #endif

         #endif
         #if KB > 187
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4376(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 188
            vbroadcastsd 4384(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 188
            vbroadcastsd 4392(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 188
            vmovapd 17920(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 188
            vmovapd 17952(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 188
            vmovapd 17984(pA), rA2
            #endif

         #endif
         #if KB > 188
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4400(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 189
            vbroadcastsd 4408(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 189
            vbroadcastsd 4416(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 189
            vmovapd 18016(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 189
            vmovapd 18048(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 189
            vmovapd 18080(pA), rA2
            #endif

         #endif
         #if KB > 189
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4424(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 190
            vbroadcastsd 4432(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 190
            vbroadcastsd 4440(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 190
            vmovapd 18112(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 190
            vmovapd 18144(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 190
            vmovapd 18176(pA), rA2
            #endif

         #endif
         #if KB > 190
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4448(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 191
            vbroadcastsd 4456(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 191
            vbroadcastsd 4464(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 191
            vmovapd 18208(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 191
            vmovapd 18240(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 191
            vmovapd 18272(pA), rA2
            #endif

         #endif
         #if KB > 191
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4472(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 192
            vbroadcastsd -128(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 192
            vbroadcastsd -120(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 192
            vmovapd 18304(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 192
            vmovapd 18336(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 192
            vmovapd 18368(pA), rA2
            #endif

         #endif
         #if KB > 192
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -112(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 193
            vbroadcastsd -104(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 193
            vbroadcastsd -96(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 193
            vmovapd 18400(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 193
            vmovapd 18432(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 193
            vmovapd 18464(pA), rA2
            #endif

         #endif
         #if KB > 193
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -88(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 194
            vbroadcastsd -80(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 194
            vbroadcastsd -72(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 194
            vmovapd 18496(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 194
            vmovapd 18528(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 194
            vmovapd 18560(pA), rA2
            #endif

         #endif
         #if KB > 194
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -64(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 195
            vbroadcastsd -56(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 195
            vbroadcastsd -48(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 195
            vmovapd 18592(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 195
            vmovapd 18624(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 195
            vmovapd 18656(pA), rA2
            #endif

         #endif
         #if KB > 195
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -40(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 196
            vbroadcastsd -32(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 196
            vbroadcastsd -24(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 196
            vmovapd 18688(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 196
            vmovapd 18720(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 196
            vmovapd 18752(pA), rA2
            #endif

         #endif
         #if KB > 196
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -16(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 197
            vbroadcastsd -8(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 197
            vbroadcastsd (pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 197
            vmovapd 18784(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 197
            vmovapd 18816(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 197
            vmovapd 18848(pA), rA2
            #endif

         #endif
         #if KB > 197
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 8(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 198
            vbroadcastsd 16(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 198
            vbroadcastsd 24(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 198
            vmovapd 18880(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 198
            vmovapd 18912(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 198
            vmovapd 18944(pA), rA2
            #endif

         #endif
         #if KB > 198
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 32(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 199
            vbroadcastsd 40(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 199
            vbroadcastsd 48(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 199
            vmovapd 18976(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 199
            vmovapd 19008(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 199
            vmovapd 19040(pA), rA2
            #endif

         #endif
         #if KB > 199
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 56(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 200
            vbroadcastsd 64(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 200
            vbroadcastsd 72(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 200
            vmovapd 19072(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 200
            vmovapd 19104(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 200
            vmovapd 19136(pA), rA2
            #endif

         #endif
         #if KB > 200
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 80(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 201
            vbroadcastsd 88(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 201
            vbroadcastsd 96(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 201
            vmovapd 19168(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 201
            vmovapd 19200(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 201
            vmovapd 19232(pA), rA2
            #endif

         #endif
         #if KB > 201
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 104(pB,i2304,2), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 202
            vbroadcastsd 112(pB,i2304,2), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 202
            vbroadcastsd 120(pB,i2304,2), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 202
            vmovapd 19264(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 202
            vmovapd 19296(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 202
            vmovapd 19328(pA), rA2
            #endif

         #endif
         #if KB > 202
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4736(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 203
            vbroadcastsd 4744(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 203
            vbroadcastsd 4752(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 203
            vmovapd 19360(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 203
            vmovapd 19392(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 203
            vmovapd 19424(pA), rA2
            #endif

         #endif
         #if KB > 203
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4760(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 204
            vbroadcastsd 4768(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 204
            vbroadcastsd 4776(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 204
            vmovapd 19456(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 204
            vmovapd 19488(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 204
            vmovapd 19520(pA), rA2
            #endif

         #endif
         #if KB > 204
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4784(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 205
            vbroadcastsd 4792(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 205
            vbroadcastsd 4800(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 205
            vmovapd 19552(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 205
            vmovapd 19584(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 205
            vmovapd 19616(pA), rA2
            #endif

         #endif
         #if KB > 205
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4808(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 206
            vbroadcastsd 4816(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 206
            vbroadcastsd 4824(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 206
            vmovapd 19648(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 206
            vmovapd 19680(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 206
            vmovapd 19712(pA), rA2
            #endif

         #endif
         #if KB > 206
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4832(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 207
            vbroadcastsd 4840(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 207
            vbroadcastsd 4848(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 207
            vmovapd 19744(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 207
            vmovapd 19776(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 207
            vmovapd 19808(pA), rA2
            #endif

         #endif
         #if KB > 207
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4856(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 208
            vbroadcastsd 4864(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 208
            vbroadcastsd 4872(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 208
            vmovapd 19840(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 208
            vmovapd 19872(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 208
            vmovapd 19904(pA), rA2
            #endif

         #endif
         #if KB > 208
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4880(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 209
            vbroadcastsd 4888(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 209
            vbroadcastsd 4896(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 209
            vmovapd 19936(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 209
            vmovapd 19968(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 209
            vmovapd 20000(pA), rA2
            #endif

         #endif
         #if KB > 209
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4904(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 210
            vbroadcastsd 4912(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 210
            vbroadcastsd 4920(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 210
            vmovapd 20032(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 210
            vmovapd 20064(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 210
            vmovapd 20096(pA), rA2
            #endif

         #endif
         #if KB > 210
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4928(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 211
            vbroadcastsd 4936(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 211
            vbroadcastsd 4944(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 211
            vmovapd 20128(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 211
            vmovapd 20160(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 211
            vmovapd 20192(pA), rA2
            #endif

         #endif
         #if KB > 211
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4952(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 212
            vbroadcastsd 4960(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 212
            vbroadcastsd 4968(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 212
            vmovapd 20224(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 212
            vmovapd 20256(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 212
            vmovapd 20288(pA), rA2
            #endif

         #endif
         #if KB > 212
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 4976(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 213
            vbroadcastsd 4984(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 213
            vbroadcastsd -128(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 213
            vmovapd 20320(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 213
            vmovapd 20352(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 213
            vmovapd 20384(pA), rA2
            #endif

         #endif
         #if KB > 213
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -120(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 214
            vbroadcastsd -112(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 214
            vbroadcastsd -104(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 214
            vmovapd 20416(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 214
            vmovapd 20448(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 214
            vmovapd 20480(pA), rA2
            #endif

         #endif
         #if KB > 214
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -96(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 215
            vbroadcastsd -88(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 215
            vbroadcastsd -80(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 215
            vmovapd 20512(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 215
            vmovapd 20544(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 215
            vmovapd 20576(pA), rA2
            #endif

         #endif
         #if KB > 215
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -72(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 216
            vbroadcastsd -64(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 216
            vbroadcastsd -56(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 216
            vmovapd 20608(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 216
            vmovapd 20640(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 216
            vmovapd 20672(pA), rA2
            #endif

         #endif
         #if KB > 216
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -48(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 217
            vbroadcastsd -40(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 217
            vbroadcastsd -32(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 217
            vmovapd 20704(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 217
            vmovapd 20736(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 217
            vmovapd 20768(pA), rA2
            #endif

         #endif
         #if KB > 217
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd -24(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 218
            vbroadcastsd -16(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 218
            vbroadcastsd -8(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 218
            vmovapd 20800(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 218
            vmovapd 20832(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 218
            vmovapd 20864(pA), rA2
            #endif

         #endif
         #if KB > 218
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd (pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 219
            vbroadcastsd 8(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 219
            vbroadcastsd 16(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 219
            vmovapd 20896(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 219
            vmovapd 20928(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 219
            vmovapd 20960(pA), rA2
            #endif

         #endif
         #if KB > 219
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 24(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 220
            vbroadcastsd 32(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 220
            vbroadcastsd 40(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 220
            vmovapd 20992(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 220
            vmovapd 21024(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 220
            vmovapd 21056(pA), rA2
            #endif

         #endif
         #if KB > 220
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 48(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 221
            vbroadcastsd 56(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 221
            vbroadcastsd 64(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 221
            vmovapd 21088(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 221
            vmovapd 21120(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 221
            vmovapd 21152(pA), rA2
            #endif

         #endif
         #if KB > 221
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 72(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 222
            vbroadcastsd 80(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 222
            vbroadcastsd 88(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 222
            vmovapd 21184(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 222
            vmovapd 21216(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 222
            vmovapd 21248(pA), rA2
            #endif

         #endif
         #if KB > 222
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 96(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 223
            vbroadcastsd 104(pB,i1280,4), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 223
            vbroadcastsd 112(pB,i1280,4), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 223
            vmovapd 21280(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 223
            vmovapd 21312(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 223
            vmovapd 21344(pA), rA2
            #endif

         #endif
         #if KB > 223
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 120(pB,i1280,4), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 224
            vbroadcastsd 5248(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 224
            vbroadcastsd 5256(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 224
            vmovapd 21376(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 224
            vmovapd 21408(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 224
            vmovapd 21440(pA), rA2
            #endif

         #endif
         #if KB > 224
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5264(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 225
            vbroadcastsd 5272(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 225
            vbroadcastsd 5280(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 225
            vmovapd 21472(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 225
            vmovapd 21504(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 225
            vmovapd 21536(pA), rA2
            #endif

         #endif
         #if KB > 225
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5288(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 226
            vbroadcastsd 5296(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 226
            vbroadcastsd 5304(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 226
            vmovapd 21568(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 226
            vmovapd 21600(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 226
            vmovapd 21632(pA), rA2
            #endif

         #endif
         #if KB > 226
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5312(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 227
            vbroadcastsd 5320(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 227
            vbroadcastsd 5328(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 227
            vmovapd 21664(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 227
            vmovapd 21696(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 227
            vmovapd 21728(pA), rA2
            #endif

         #endif
         #if KB > 227
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5336(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 228
            vbroadcastsd 5344(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 228
            vbroadcastsd 5352(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 228
            vmovapd 21760(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 228
            vmovapd 21792(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 228
            vmovapd 21824(pA), rA2
            #endif

         #endif
         #if KB > 228
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5360(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 229
            vbroadcastsd 5368(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 229
            vbroadcastsd 5376(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 229
            vmovapd 21856(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 229
            vmovapd 21888(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 229
            vmovapd 21920(pA), rA2
            #endif

         #endif
         #if KB > 229
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5384(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 230
            vbroadcastsd 5392(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 230
            vbroadcastsd 5400(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 230
            vmovapd 21952(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 230
            vmovapd 21984(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 230
            vmovapd 22016(pA), rA2
            #endif

         #endif
         #if KB > 230
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5408(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 231
            vbroadcastsd 5416(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 231
            vbroadcastsd 5424(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 231
            vmovapd 22048(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 231
            vmovapd 22080(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 231
            vmovapd 22112(pA), rA2
            #endif

         #endif
         #if KB > 231
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5432(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 232
            vbroadcastsd 5440(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 232
            vbroadcastsd 5448(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 232
            vmovapd 22144(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 232
            vmovapd 22176(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 232
            vmovapd 22208(pA), rA2
            #endif

         #endif
         #if KB > 232
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5456(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 233
            vbroadcastsd 5464(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 233
            vbroadcastsd 5472(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 233
            vmovapd 22240(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 233
            vmovapd 22272(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 233
            vmovapd 22304(pA), rA2
            #endif

         #endif
         #if KB > 233
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5480(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 234
            vbroadcastsd 5488(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 234
            vbroadcastsd 5496(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 234
            vmovapd 22336(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 234
            vmovapd 22368(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 234
            vmovapd 22400(pA), rA2
            #endif

         #endif
         #if KB > 234
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5504(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 235
            vbroadcastsd 5512(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 235
            vbroadcastsd 5520(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 235
            vmovapd 22432(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 235
            vmovapd 22464(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 235
            vmovapd 22496(pA), rA2
            #endif

         #endif
         #if KB > 235
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5528(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 236
            vbroadcastsd 5536(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 236
            vbroadcastsd 5544(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 236
            vmovapd 22528(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 236
            vmovapd 22560(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 236
            vmovapd 22592(pA), rA2
            #endif

         #endif
         #if KB > 236
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5552(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 237
            vbroadcastsd 5560(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 237
            vbroadcastsd 5568(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 237
            vmovapd 22624(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 237
            vmovapd 22656(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 237
            vmovapd 22688(pA), rA2
            #endif

         #endif
         #if KB > 237
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5576(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 238
            vbroadcastsd 5584(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 238
            vbroadcastsd 5592(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 238
            vmovapd 22720(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 238
            vmovapd 22752(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 238
            vmovapd 22784(pA), rA2
            #endif

         #endif
         #if KB > 238
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5600(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 239
            vbroadcastsd 5608(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 239
            vbroadcastsd 5616(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 239
            vmovapd 22816(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 239
            vmovapd 22848(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 239
            vmovapd 22880(pA), rA2
            #endif

         #endif
         #if KB > 239
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5624(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 240
            vbroadcastsd 5632(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 240
            vbroadcastsd 5640(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 240
            vmovapd 22912(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 240
            vmovapd 22944(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 240
            vmovapd 22976(pA), rA2
            #endif

         #endif
         #if KB > 240
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5648(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 241
            vbroadcastsd 5656(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 241
            vbroadcastsd 5664(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 241
            vmovapd 23008(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 241
            vmovapd 23040(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 241
            vmovapd 23072(pA), rA2
            #endif

         #endif
         #if KB > 241
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5672(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 242
            vbroadcastsd 5680(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 242
            vbroadcastsd 5688(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 242
            vmovapd 23104(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 242
            vmovapd 23136(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 242
            vmovapd 23168(pA), rA2
            #endif

         #endif
         #if KB > 242
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5696(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 243
            vbroadcastsd 5704(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 243
            vbroadcastsd 5712(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 243
            vmovapd 23200(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 243
            vmovapd 23232(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 243
            vmovapd 23264(pA), rA2
            #endif

         #endif
         #if KB > 243
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5720(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 244
            vbroadcastsd 5728(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 244
            vbroadcastsd 5736(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 244
            vmovapd 23296(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 244
            vmovapd 23328(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 244
            vmovapd 23360(pA), rA2
            #endif

         #endif
         #if KB > 244
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5744(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 245
            vbroadcastsd 5752(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 245
            vbroadcastsd 5760(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 245
            vmovapd 23392(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 245
            vmovapd 23424(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 245
            vmovapd 23456(pA), rA2
            #endif

         #endif
         #if KB > 245
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5768(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 246
            vbroadcastsd 5776(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 246
            vbroadcastsd 5784(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 246
            vmovapd 23488(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 246
            vmovapd 23520(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 246
            vmovapd 23552(pA), rA2
            #endif

         #endif
         #if KB > 246
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5792(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 247
            vbroadcastsd 5800(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 247
            vbroadcastsd 5808(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 247
            vmovapd 23584(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 247
            vmovapd 23616(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 247
            vmovapd 23648(pA), rA2
            #endif

         #endif
         #if KB > 247
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5816(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 248
            vbroadcastsd 5824(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 248
            vbroadcastsd 5832(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 248
            vmovapd 23680(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 248
            vmovapd 23712(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 248
            vmovapd 23744(pA), rA2
            #endif

         #endif
         #if KB > 248
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5840(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 249
            vbroadcastsd 5848(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 249
            vbroadcastsd 5856(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 249
            vmovapd 23776(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 249
            vmovapd 23808(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 249
            vmovapd 23840(pA), rA2
            #endif

         #endif
         #if KB > 249
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5864(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 250
            vbroadcastsd 5872(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 250
            vbroadcastsd 5880(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 250
            vmovapd 23872(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 250
            vmovapd 23904(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 250
            vmovapd 23936(pA), rA2
            #endif

         #endif
         #if KB > 250
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5888(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 251
            vbroadcastsd 5896(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 251
            vbroadcastsd 5904(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 251
            vmovapd 23968(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 251
            vmovapd 24000(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 251
            vmovapd 24032(pA), rA2
            #endif

         #endif
         #if KB > 251
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5912(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 252
            vbroadcastsd 5920(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 252
            vbroadcastsd 5928(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 252
            vmovapd 24064(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 252
            vmovapd 24096(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 252
            vmovapd 24128(pA), rA2
            #endif

         #endif
         #if KB > 252
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5936(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 253
            vbroadcastsd 5944(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 253
            vbroadcastsd 5952(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 253
            vmovapd 24160(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 253
            vmovapd 24192(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 253
            vmovapd 24224(pA), rA2
            #endif

         #endif
         #if KB > 253
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5960(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 254
            vbroadcastsd 5968(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 254
            vbroadcastsd 5976(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 254
            vmovapd 24256(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 254
            vmovapd 24288(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 254
            vmovapd 24320(pA), rA2
            #endif

         #endif
         #if KB > 254
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 5984(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 255
            vbroadcastsd 5992(pB), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 255
            vbroadcastsd 6000(pB), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 255
            vmovapd 24352(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 255
            vmovapd 24384(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 255
            vmovapd 24416(pA), rA2
            #endif

         #endif
         #if KB > 255
            vmulpd rB0, rA0, m0
            vaddpd rC00, m0, rC00
            vbroadcastsd 6008(pB), rB2
            vmulpd rB0, rA1, m0
            vaddpd rC10, m0, rC10
            vmulpd rB0, rA2, m0
            vaddpd rC20, m0, rC20
            #if KB > 256
            vbroadcastsd -128(pB,i768,8), rB0
            #endif

            vmulpd rB1, rA0, m0
            vaddpd rC01, m0, rC01
            vmulpd rB1, rA1, m0
            vaddpd rC11, m0, rC11
            vmulpd rB1, rA2, m0
            vaddpd rC21, m0, rC21
            #if KB > 256
            vbroadcastsd -120(pB,i768,8), rB1
            #endif

            vmulpd rB2, rA0, m0
            vaddpd rC02, m0, rC02
            #if KB > 256
            vmovapd 24448(pA), rA0
            #endif
            vmulpd rB2, rA1, m0
            vaddpd rC12, m0, rC12
            #if KB > 256
            vmovapd 24480(pA), rA1
            #endif
            vmulpd rB2, rA2, m0
            vaddpd rC22, m0, rC22
            #if KB > 256
            vmovapd 24512(pA), rA2
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
         add $KB*3*8, pB        /* pB += K*NU*sizeof */
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
