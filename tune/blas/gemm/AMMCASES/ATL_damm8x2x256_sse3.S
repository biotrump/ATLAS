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
#elif defined(BETA1)
   #define BETCOP addpd
#endif
#define MOVAPD movaps
#define UNPCKHPD movhlps

#define m0      %xmm0
#define rA0     %xmm1
#define rA1     %xmm2
#define rA2     %xmm3
#define rA3     %xmm4
#define rB0     %xmm5
#define rB1     %xmm6
#define rb1     %xmm7
#define rC00    %xmm8
#define rC10    %xmm9
#define rC20    %xmm10
#define rC30    %xmm11
#define rC01    %xmm12
#define rC11    %xmm13
#define rC21    %xmm14
#define rC31    %xmm15

#define pA      %rcx
#define pB      %rdx
#define pC      %rbp
#define pfA     %rax
#define pfB     %rbx
#define pfC     %rdi
#define nnu     %rsi

#define incPF   %r8
#define incA    %r9
#define incB    %r10
#define pA0     %r11
#define nnu0    %r12
#define r256    %r13
#define nmu     %r14
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
   #define FSIZE 6*8
/*
 * Save callee-saved iregs
 */
   sub $FSIZE, %rsp
   movq    %rbp, (%rsp)
   movq    %rbx, 8(%rsp)
   movq    %r12, 16(%rsp)
   movq    %r13, 24(%rsp)
   movq    %r14, 32(%rsp)
/*
 * Load parameters
 */
   mov  %rdi, nmu
   mov  %rsi, nnu0
   mov  %r8, pB
   mov  %r9, pC
   movq FSIZE+8(%rsp), pfA
   movq FSIZE+16(%rsp), pfB
   movq FSIZE+24(%rsp), pfC
   mov  $256, r256
   mov  $8*8*2, incPF
/*
 * Add 128 to ptrs to maximize unrolling range
 */
   add incPF, pA
   add incPF, pB
   mov pA, pA0
   mov pB, pB0
   mov nnu0, nnu
   mov $2*KB*8, incB   /* incB = nu*KB*sizeof */
   mov $8*KB*8, incA   /* incA = mu*KB*sizeof */

#if 0
   MLOOP:
      mov nnu0, nnu
      NLOOP:
#else
   ALIGN16
   .local MNLOOP
   MNLOOP:
#endif
/*
 *       Unroll first iteration to avoid zeroing rCxx
 */
         MOVAPD -128(pB), rB1
         movddup rB1, rB0
         MOVAPD -128(pA), rC00
         MOVAPD rC00, rC01
         mulpd rB0, rC00
         MOVAPD -112(pA), rC10
         MOVAPD rC10, rC11
         mulpd rB0, rC10
         MOVAPD -96(pA), rC20
         MOVAPD rC20, rC21
         mulpd rB0, rC20
         MOVAPD -80(pA), rC30
         MOVAPD rC30, rC31
         mulpd rB0, rC30


         UNPCKHPD rB1, rB1
         #if KB > 1
            MOVAPD -112(pB), rb1
         #endif
         mulpd rB1, rC01
         #if KB > 1
            movddup rb1, rB0
            MOVAPD -64(pA), rA0
         #endif
         mulpd rB1, rC11
         #if KB > 1
            MOVAPD -48(pA), rA1
         #endif
         mulpd rB1, rC21
         #if KB > 1
            MOVAPD -32(pA), rA2
         #endif
         mulpd rB1, rC31
         #if KB > 1
            UNPCKHPD rb1, rb1
         #endif

         #if KB > 1
               MOVAPD -16(pA), rA3
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 2
               MOVAPD -96(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 2
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 2
               MOVAPD (pA), rA0
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 2
               MOVAPD 16(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 2
               MOVAPD 32(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 2
               MOVAPD 48(pA), rA3
               UNPCKHPD rB1, rB1
            #endif
         #endif

         #if KB > 2
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 3
               MOVAPD -80(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 3
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 3
               UNPCKHPD rb1, rb1
               MOVAPD 64(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 3
               MOVAPD 80(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 3
               MOVAPD 96(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 3
               MOVAPD 112(pA), rA3
              add r256, pA
            #endif
         #endif

         #if KB > 3
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 4
               MOVAPD -64(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
               prefC 64(pC)
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
               prefC 64(pC)
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 4
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 4
               MOVAPD -128(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 4
               MOVAPD -112(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 4
               MOVAPD -96(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 4
               MOVAPD -80(pA), rA3
            #endif
         #endif

         ALIGN16
         #if KB > 4
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 5
               MOVAPD -48(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 5
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 5
               UNPCKHPD rb1, rb1
               MOVAPD -64(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 5
               MOVAPD -48(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 5
               MOVAPD -32(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 5
               MOVAPD -16(pA), rA3
            #endif
         #endif

         #if KB > 5
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 6
               MOVAPD -32(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 6
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 6
               MOVAPD 0(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 6
               MOVAPD 16(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 6
               MOVAPD 32(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 6
               MOVAPD 48(pA), rA3
            #endif
         #endif
         #if KB > 6
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 7
               MOVAPD -16(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 7
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 7
               UNPCKHPD rb1, rb1
               MOVAPD 64(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 7
               MOVAPD 80(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 7
               MOVAPD 96(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 7
               MOVAPD 112(pA), rA3
            #endif
         #endif

         #if KB > 7
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 8
               MOVAPD 0(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 8
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 8
               MOVAPD 128(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 8
               MOVAPD 144(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 8
               MOVAPD 160(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 8
               MOVAPD 176(pA), rA3
            #endif
         #endif
         #if KB > 8
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 9
               MOVAPD 16(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 9
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 9
               UNPCKHPD rb1, rb1
               MOVAPD 192(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 9
               MOVAPD 208(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 9
               MOVAPD 224(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 9
               MOVAPD 240(pA), rA3
            #endif
         #endif

         #if KB > 9
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 10
               MOVAPD 32(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 10
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 10
               MOVAPD 256(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 10
               MOVAPD 272(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 10
               MOVAPD 288(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 10
               MOVAPD 304(pA), rA3
            #endif
         #endif
         #if KB > 10
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 11
               MOVAPD 48(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 11
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 11
               UNPCKHPD rb1, rb1
               MOVAPD 320(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 11
               MOVAPD 336(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 11
               MOVAPD 352(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 11
               MOVAPD 368(pA), rA3
            #endif
         #endif

         #if KB > 11
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 12
               MOVAPD 64(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 12
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 12
               MOVAPD 384(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 12
               MOVAPD 400(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 12
               MOVAPD 416(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 12
               MOVAPD 432(pA), rA3
            #endif
         #endif
         #if KB > 12
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 13
               MOVAPD 80(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 13
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 13
               UNPCKHPD rb1, rb1
               MOVAPD 448(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 13
               MOVAPD 464(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 13
               MOVAPD 480(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 13
               MOVAPD 496(pA), rA3
            #endif
         #endif

         #if KB > 13
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 14
               MOVAPD 96(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 14
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 14
               MOVAPD 512(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 14
               MOVAPD 528(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 14
               MOVAPD 544(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 14
               MOVAPD 560(pA), rA3
            #endif
         #endif
         #if KB > 14
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 15
               MOVAPD 112(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 15
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 15
               UNPCKHPD rb1, rb1
               MOVAPD 576(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 15
               MOVAPD 592(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 15
               MOVAPD 608(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 15
               MOVAPD 624(pA), rA3
            #endif
         #endif

         #if KB > 15
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 16
               MOVAPD 128(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 16
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 16
               MOVAPD 640(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 16
               MOVAPD 656(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 16
               MOVAPD 672(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 16
               MOVAPD 688(pA), rA3
            #endif
         #endif
         #if KB > 16
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 17
               MOVAPD 144(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 17
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 17
               UNPCKHPD rb1, rb1
               MOVAPD 704(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 17
               MOVAPD 720(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 17
               MOVAPD 736(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 17
               MOVAPD 752(pA), rA3
            #endif
         #endif

         #if KB > 17
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 18
               MOVAPD 160(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 18
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 18
               MOVAPD 768(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 18
               MOVAPD 784(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 18
               MOVAPD 800(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 18
               MOVAPD 816(pA), rA3
            #endif
         #endif
         #if KB > 18
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 19
               MOVAPD 176(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 19
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 19
               UNPCKHPD rb1, rb1
               MOVAPD 832(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 19
               MOVAPD 848(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 19
               MOVAPD 864(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 19
               MOVAPD 880(pA), rA3
            #endif
         #endif

         #if KB > 19
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 20
               MOVAPD 192(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 20
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 20
               MOVAPD 896(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 20
               MOVAPD 912(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 20
               MOVAPD 928(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 20
               MOVAPD 944(pA), rA3
            #endif
         #endif
         #if KB > 20
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 21
               MOVAPD 208(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 21
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 21
               UNPCKHPD rb1, rb1
               MOVAPD 960(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 21
               MOVAPD 976(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 21
               MOVAPD 992(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 21
               MOVAPD 1008(pA), rA3
            #endif
         #endif

         #if KB > 21
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 22
               MOVAPD 224(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 22
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 22
               MOVAPD 1024(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 22
               MOVAPD 1040(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 22
               MOVAPD 1056(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 22
               MOVAPD 1072(pA), rA3
            #endif
         #endif
         #if KB > 22
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 23
               MOVAPD 240(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 23
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 23
               UNPCKHPD rb1, rb1
               MOVAPD 1088(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 23
               MOVAPD 1104(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 23
               MOVAPD 1120(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 23
               MOVAPD 1136(pA), rA3
            #endif
         #endif

         #if KB > 23
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 24
               MOVAPD 256(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 24
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 24
               MOVAPD 1152(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 24
               MOVAPD 1168(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 24
               MOVAPD 1184(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 24
               MOVAPD 1200(pA), rA3
            #endif
         #endif
         #if KB > 24
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 25
               MOVAPD 272(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 25
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 25
               UNPCKHPD rb1, rb1
               MOVAPD 1216(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 25
               MOVAPD 1232(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 25
               MOVAPD 1248(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 25
               MOVAPD 1264(pA), rA3
            #endif
         #endif

         #if KB > 25
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 26
               MOVAPD 288(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 26
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 26
               MOVAPD 1280(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 26
               MOVAPD 1296(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 26
               MOVAPD 1312(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 26
               MOVAPD 1328(pA), rA3
            #endif
         #endif
         #if KB > 26
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 27
               MOVAPD 304(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 27
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 27
               UNPCKHPD rb1, rb1
               MOVAPD 1344(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 27
               MOVAPD 1360(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 27
               MOVAPD 1376(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 27
               MOVAPD 1392(pA), rA3
            #endif
         #endif

         #if KB > 27
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 28
               MOVAPD 320(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 28
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 28
               MOVAPD 1408(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 28
               MOVAPD 1424(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 28
               MOVAPD 1440(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 28
               MOVAPD 1456(pA), rA3
            #endif
         #endif
         #if KB > 28
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 29
               MOVAPD 336(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 29
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 29
               UNPCKHPD rb1, rb1
               MOVAPD 1472(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 29
               MOVAPD 1488(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 29
               MOVAPD 1504(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 29
               MOVAPD 1520(pA), rA3
            #endif
         #endif

         #if KB > 29
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 30
               MOVAPD 352(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 30
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 30
               MOVAPD 1536(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 30
               MOVAPD 1552(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 30
               MOVAPD 1568(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 30
               MOVAPD 1584(pA), rA3
            #endif
         #endif
         #if KB > 30
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 31
               MOVAPD 368(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 31
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 31
               UNPCKHPD rb1, rb1
               MOVAPD 1600(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 31
               MOVAPD 1616(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 31
               MOVAPD 1632(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 31
               MOVAPD 1648(pA), rA3
            #endif
         #endif

         #if KB > 31
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 32
               MOVAPD 384(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 32
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 32
               MOVAPD 1664(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 32
               MOVAPD 1680(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 32
               MOVAPD 1696(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 32
               MOVAPD 1712(pA), rA3
            #endif
         #endif
         #if KB > 32
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 33
               MOVAPD 400(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 33
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 33
               UNPCKHPD rb1, rb1
               MOVAPD 1728(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 33
               MOVAPD 1744(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 33
               MOVAPD 1760(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 33
               MOVAPD 1776(pA), rA3
            #endif
         #endif

         #if KB > 33
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 34
               MOVAPD 416(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 34
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 34
               MOVAPD 1792(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 34
               MOVAPD 1808(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 34
               MOVAPD 1824(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 34
               MOVAPD 1840(pA), rA3
            #endif
         #endif
         #if KB > 34
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 35
               MOVAPD 432(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 35
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 35
               UNPCKHPD rb1, rb1
               MOVAPD 1856(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 35
               MOVAPD 1872(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 35
               MOVAPD 1888(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 35
               MOVAPD 1904(pA), rA3
            #endif
         #endif

         #if KB > 35
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 36
               MOVAPD 448(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 36
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 36
               MOVAPD 1920(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 36
               MOVAPD 1936(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 36
               MOVAPD 1952(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 36
               MOVAPD 1968(pA), rA3
            #endif
         #endif
         #if KB > 36
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 37
               MOVAPD 464(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 37
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 37
               UNPCKHPD rb1, rb1
               MOVAPD 1984(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 37
               MOVAPD 2000(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 37
               MOVAPD 2016(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 37
               MOVAPD 2032(pA), rA3
            #endif
         #endif

         #if KB > 37
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 38
               MOVAPD 480(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 38
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 38
               MOVAPD 2048(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 38
               MOVAPD 2064(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 38
               MOVAPD 2080(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 38
               MOVAPD 2096(pA), rA3
            #endif
         #endif
         #if KB > 38
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 39
               MOVAPD 496(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 39
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 39
               UNPCKHPD rb1, rb1
               MOVAPD 2112(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 39
               MOVAPD 2128(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 39
               MOVAPD 2144(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 39
               MOVAPD 2160(pA), rA3
            #endif
         #endif

         #if KB > 39
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 40
               MOVAPD 512(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 40
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 40
               MOVAPD 2176(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 40
               MOVAPD 2192(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 40
               MOVAPD 2208(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 40
               MOVAPD 2224(pA), rA3
            #endif
         #endif
         #if KB > 40
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 41
               MOVAPD 528(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 41
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 41
               UNPCKHPD rb1, rb1
               MOVAPD 2240(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 41
               MOVAPD 2256(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 41
               MOVAPD 2272(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 41
               MOVAPD 2288(pA), rA3
            #endif
         #endif

         #if KB > 41
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 42
               MOVAPD 544(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 42
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 42
               MOVAPD 2304(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 42
               MOVAPD 2320(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 42
               MOVAPD 2336(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 42
               MOVAPD 2352(pA), rA3
            #endif
         #endif
         #if KB > 42
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 43
               MOVAPD 560(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 43
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 43
               UNPCKHPD rb1, rb1
               MOVAPD 2368(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 43
               MOVAPD 2384(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 43
               MOVAPD 2400(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 43
               MOVAPD 2416(pA), rA3
            #endif
         #endif

         #if KB > 43
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 44
               MOVAPD 576(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 44
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 44
               MOVAPD 2432(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 44
               MOVAPD 2448(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 44
               MOVAPD 2464(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 44
               MOVAPD 2480(pA), rA3
            #endif
         #endif
         #if KB > 44
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 45
               MOVAPD 592(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 45
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 45
               UNPCKHPD rb1, rb1
               MOVAPD 2496(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 45
               MOVAPD 2512(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 45
               MOVAPD 2528(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 45
               MOVAPD 2544(pA), rA3
            #endif
         #endif

         #if KB > 45
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 46
               MOVAPD 608(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 46
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 46
               MOVAPD 2560(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 46
               MOVAPD 2576(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 46
               MOVAPD 2592(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 46
               MOVAPD 2608(pA), rA3
            #endif
         #endif
         #if KB > 46
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 47
               MOVAPD 624(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 47
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 47
               UNPCKHPD rb1, rb1
               MOVAPD 2624(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 47
               MOVAPD 2640(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 47
               MOVAPD 2656(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 47
               MOVAPD 2672(pA), rA3
            #endif
         #endif

         #if KB > 47
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 48
               MOVAPD 640(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 48
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 48
               MOVAPD 2688(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 48
               MOVAPD 2704(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 48
               MOVAPD 2720(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 48
               MOVAPD 2736(pA), rA3
            #endif
         #endif
         #if KB > 48
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 49
               MOVAPD 656(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 49
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 49
               UNPCKHPD rb1, rb1
               MOVAPD 2752(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 49
               MOVAPD 2768(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 49
               MOVAPD 2784(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 49
               MOVAPD 2800(pA), rA3
            #endif
         #endif

         #if KB > 49
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 50
               MOVAPD 672(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 50
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 50
               MOVAPD 2816(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 50
               MOVAPD 2832(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 50
               MOVAPD 2848(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 50
               MOVAPD 2864(pA), rA3
            #endif
         #endif
         #if KB > 50
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 51
               MOVAPD 688(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 51
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 51
               UNPCKHPD rb1, rb1
               MOVAPD 2880(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 51
               MOVAPD 2896(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 51
               MOVAPD 2912(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 51
               MOVAPD 2928(pA), rA3
            #endif
         #endif

         #if KB > 51
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 52
               MOVAPD 704(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 52
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 52
               MOVAPD 2944(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 52
               MOVAPD 2960(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 52
               MOVAPD 2976(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 52
               MOVAPD 2992(pA), rA3
            #endif
         #endif
         #if KB > 52
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 53
               MOVAPD 720(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 53
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 53
               UNPCKHPD rb1, rb1
               MOVAPD 3008(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 53
               MOVAPD 3024(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 53
               MOVAPD 3040(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 53
               MOVAPD 3056(pA), rA3
            #endif
         #endif

         #if KB > 53
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 54
               MOVAPD 736(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 54
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 54
               MOVAPD 3072(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 54
               MOVAPD 3088(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 54
               MOVAPD 3104(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 54
               MOVAPD 3120(pA), rA3
            #endif
         #endif
         #if KB > 54
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 55
               MOVAPD 752(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 55
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 55
               UNPCKHPD rb1, rb1
               MOVAPD 3136(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 55
               MOVAPD 3152(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 55
               MOVAPD 3168(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 55
               MOVAPD 3184(pA), rA3
            #endif
         #endif

         #if KB > 55
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 56
               MOVAPD 768(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 56
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 56
               MOVAPD 3200(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 56
               MOVAPD 3216(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 56
               MOVAPD 3232(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 56
               MOVAPD 3248(pA), rA3
            #endif
         #endif
         #if KB > 56
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 57
               MOVAPD 784(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 57
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 57
               UNPCKHPD rb1, rb1
               MOVAPD 3264(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 57
               MOVAPD 3280(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 57
               MOVAPD 3296(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 57
               MOVAPD 3312(pA), rA3
            #endif
         #endif

         #if KB > 57
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 58
               MOVAPD 800(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 58
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 58
               MOVAPD 3328(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 58
               MOVAPD 3344(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 58
               MOVAPD 3360(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 58
               MOVAPD 3376(pA), rA3
            #endif
         #endif
         #if KB > 58
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 59
               MOVAPD 816(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 59
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 59
               UNPCKHPD rb1, rb1
               MOVAPD 3392(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 59
               MOVAPD 3408(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 59
               MOVAPD 3424(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 59
               MOVAPD 3440(pA), rA3
            #endif
         #endif

         #if KB > 59
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 60
               MOVAPD 832(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 60
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 60
               MOVAPD 3456(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 60
               MOVAPD 3472(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 60
               MOVAPD 3488(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 60
               MOVAPD 3504(pA), rA3
            #endif
         #endif
         #if KB > 60
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 61
               MOVAPD 848(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 61
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 61
               UNPCKHPD rb1, rb1
               MOVAPD 3520(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 61
               MOVAPD 3536(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 61
               MOVAPD 3552(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 61
               MOVAPD 3568(pA), rA3
            #endif
         #endif

         #if KB > 61
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 62
               MOVAPD 864(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 62
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 62
               MOVAPD 3584(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 62
               MOVAPD 3600(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 62
               MOVAPD 3616(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 62
               MOVAPD 3632(pA), rA3
            #endif
         #endif
         #if KB > 62
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 63
               MOVAPD 880(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 63
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 63
               UNPCKHPD rb1, rb1
               MOVAPD 3648(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 63
               MOVAPD 3664(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 63
               MOVAPD 3680(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 63
               MOVAPD 3696(pA), rA3
            #endif
         #endif

         #if KB > 63
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 64
               MOVAPD 896(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 64
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 64
               MOVAPD 3712(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 64
               MOVAPD 3728(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 64
               MOVAPD 3744(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 64
               MOVAPD 3760(pA), rA3
            #endif
         #endif
         #if KB > 64
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 65
               MOVAPD 912(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 65
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 65
               UNPCKHPD rb1, rb1
               MOVAPD 3776(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 65
               MOVAPD 3792(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 65
               MOVAPD 3808(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 65
               MOVAPD 3824(pA), rA3
            #endif
         #endif

         #if KB > 65
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 66
               MOVAPD 928(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 66
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 66
               MOVAPD 3840(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 66
               MOVAPD 3856(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 66
               MOVAPD 3872(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 66
               MOVAPD 3888(pA), rA3
            #endif
         #endif
         #if KB > 66
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 67
               MOVAPD 944(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 67
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 67
               UNPCKHPD rb1, rb1
               MOVAPD 3904(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 67
               MOVAPD 3920(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 67
               MOVAPD 3936(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 67
               MOVAPD 3952(pA), rA3
            #endif
         #endif

         #if KB > 67
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 68
               MOVAPD 960(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 68
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 68
               MOVAPD 3968(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 68
               MOVAPD 3984(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 68
               MOVAPD 4000(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 68
               MOVAPD 4016(pA), rA3
            #endif
         #endif
         #if KB > 68
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 69
               MOVAPD 976(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 69
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 69
               UNPCKHPD rb1, rb1
               MOVAPD 4032(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 69
               MOVAPD 4048(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 69
               MOVAPD 4064(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 69
               MOVAPD 4080(pA), rA3
            #endif
         #endif

         #if KB > 69
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 70
               MOVAPD 992(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 70
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 70
               MOVAPD 4096(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 70
               MOVAPD 4112(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 70
               MOVAPD 4128(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 70
               MOVAPD 4144(pA), rA3
            #endif
         #endif
         #if KB > 70
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 71
               MOVAPD 1008(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 71
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 71
               UNPCKHPD rb1, rb1
               MOVAPD 4160(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 71
               MOVAPD 4176(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 71
               MOVAPD 4192(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 71
               MOVAPD 4208(pA), rA3
            #endif
         #endif

         #if KB > 71
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 72
               MOVAPD 1024(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 72
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 72
               MOVAPD 4224(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 72
               MOVAPD 4240(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 72
               MOVAPD 4256(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 72
               MOVAPD 4272(pA), rA3
            #endif
         #endif
         #if KB > 72
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 73
               MOVAPD 1040(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 73
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 73
               UNPCKHPD rb1, rb1
               MOVAPD 4288(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 73
               MOVAPD 4304(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 73
               MOVAPD 4320(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 73
               MOVAPD 4336(pA), rA3
            #endif
         #endif

         #if KB > 73
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 74
               MOVAPD 1056(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 74
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 74
               MOVAPD 4352(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 74
               MOVAPD 4368(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 74
               MOVAPD 4384(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 74
               MOVAPD 4400(pA), rA3
            #endif
         #endif
         #if KB > 74
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 75
               MOVAPD 1072(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 75
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 75
               UNPCKHPD rb1, rb1
               MOVAPD 4416(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 75
               MOVAPD 4432(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 75
               MOVAPD 4448(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 75
               MOVAPD 4464(pA), rA3
            #endif
         #endif

         #if KB > 75
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 76
               MOVAPD 1088(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 76
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 76
               MOVAPD 4480(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 76
               MOVAPD 4496(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 76
               MOVAPD 4512(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 76
               MOVAPD 4528(pA), rA3
            #endif
         #endif
         #if KB > 76
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 77
               MOVAPD 1104(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 77
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 77
               UNPCKHPD rb1, rb1
               MOVAPD 4544(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 77
               MOVAPD 4560(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 77
               MOVAPD 4576(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 77
               MOVAPD 4592(pA), rA3
            #endif
         #endif

         #if KB > 77
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 78
               MOVAPD 1120(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 78
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 78
               MOVAPD 4608(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 78
               MOVAPD 4624(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 78
               MOVAPD 4640(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 78
               MOVAPD 4656(pA), rA3
            #endif
         #endif
         #if KB > 78
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 79
               MOVAPD 1136(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 79
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 79
               UNPCKHPD rb1, rb1
               MOVAPD 4672(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 79
               MOVAPD 4688(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 79
               MOVAPD 4704(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 79
               MOVAPD 4720(pA), rA3
            #endif
         #endif

         #if KB > 79
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 80
               MOVAPD 1152(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 80
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 80
               MOVAPD 4736(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 80
               MOVAPD 4752(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 80
               MOVAPD 4768(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 80
               MOVAPD 4784(pA), rA3
            #endif
         #endif
         #if KB > 80
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 81
               MOVAPD 1168(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 81
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 81
               UNPCKHPD rb1, rb1
               MOVAPD 4800(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 81
               MOVAPD 4816(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 81
               MOVAPD 4832(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 81
               MOVAPD 4848(pA), rA3
            #endif
         #endif

         #if KB > 81
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 82
               MOVAPD 1184(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 82
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 82
               MOVAPD 4864(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 82
               MOVAPD 4880(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 82
               MOVAPD 4896(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 82
               MOVAPD 4912(pA), rA3
            #endif
         #endif
         #if KB > 82
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 83
               MOVAPD 1200(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 83
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 83
               UNPCKHPD rb1, rb1
               MOVAPD 4928(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 83
               MOVAPD 4944(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 83
               MOVAPD 4960(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 83
               MOVAPD 4976(pA), rA3
            #endif
         #endif

         #if KB > 83
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 84
               MOVAPD 1216(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 84
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 84
               MOVAPD 4992(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 84
               MOVAPD 5008(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 84
               MOVAPD 5024(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 84
               MOVAPD 5040(pA), rA3
            #endif
         #endif
         #if KB > 84
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 85
               MOVAPD 1232(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 85
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 85
               UNPCKHPD rb1, rb1
               MOVAPD 5056(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 85
               MOVAPD 5072(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 85
               MOVAPD 5088(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 85
               MOVAPD 5104(pA), rA3
            #endif
         #endif

         #if KB > 85
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 86
               MOVAPD 1248(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 86
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 86
               MOVAPD 5120(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 86
               MOVAPD 5136(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 86
               MOVAPD 5152(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 86
               MOVAPD 5168(pA), rA3
            #endif
         #endif
         #if KB > 86
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 87
               MOVAPD 1264(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 87
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 87
               UNPCKHPD rb1, rb1
               MOVAPD 5184(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 87
               MOVAPD 5200(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 87
               MOVAPD 5216(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 87
               MOVAPD 5232(pA), rA3
            #endif
         #endif

         #if KB > 87
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 88
               MOVAPD 1280(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 88
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 88
               MOVAPD 5248(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 88
               MOVAPD 5264(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 88
               MOVAPD 5280(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 88
               MOVAPD 5296(pA), rA3
            #endif
         #endif
         #if KB > 88
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 89
               MOVAPD 1296(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 89
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 89
               UNPCKHPD rb1, rb1
               MOVAPD 5312(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 89
               MOVAPD 5328(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 89
               MOVAPD 5344(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 89
               MOVAPD 5360(pA), rA3
            #endif
         #endif

         #if KB > 89
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 90
               MOVAPD 1312(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 90
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 90
               MOVAPD 5376(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 90
               MOVAPD 5392(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 90
               MOVAPD 5408(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 90
               MOVAPD 5424(pA), rA3
            #endif
         #endif
         #if KB > 90
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 91
               MOVAPD 1328(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 91
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 91
               UNPCKHPD rb1, rb1
               MOVAPD 5440(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 91
               MOVAPD 5456(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 91
               MOVAPD 5472(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 91
               MOVAPD 5488(pA), rA3
            #endif
         #endif

         #if KB > 91
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 92
               MOVAPD 1344(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 92
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 92
               MOVAPD 5504(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 92
               MOVAPD 5520(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 92
               MOVAPD 5536(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 92
               MOVAPD 5552(pA), rA3
            #endif
         #endif
         #if KB > 92
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 93
               MOVAPD 1360(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 93
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 93
               UNPCKHPD rb1, rb1
               MOVAPD 5568(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 93
               MOVAPD 5584(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 93
               MOVAPD 5600(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 93
               MOVAPD 5616(pA), rA3
            #endif
         #endif

         #if KB > 93
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 94
               MOVAPD 1376(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 94
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 94
               MOVAPD 5632(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 94
               MOVAPD 5648(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 94
               MOVAPD 5664(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 94
               MOVAPD 5680(pA), rA3
            #endif
         #endif
         #if KB > 94
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 95
               MOVAPD 1392(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 95
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 95
               UNPCKHPD rb1, rb1
               MOVAPD 5696(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 95
               MOVAPD 5712(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 95
               MOVAPD 5728(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 95
               MOVAPD 5744(pA), rA3
            #endif
         #endif

         #if KB > 95
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 96
               MOVAPD 1408(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 96
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 96
               MOVAPD 5760(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 96
               MOVAPD 5776(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 96
               MOVAPD 5792(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 96
               MOVAPD 5808(pA), rA3
            #endif
         #endif
         #if KB > 96
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 97
               MOVAPD 1424(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 97
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 97
               UNPCKHPD rb1, rb1
               MOVAPD 5824(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 97
               MOVAPD 5840(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 97
               MOVAPD 5856(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 97
               MOVAPD 5872(pA), rA3
            #endif
         #endif

         #if KB > 97
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 98
               MOVAPD 1440(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 98
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 98
               MOVAPD 5888(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 98
               MOVAPD 5904(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 98
               MOVAPD 5920(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 98
               MOVAPD 5936(pA), rA3
            #endif
         #endif
         #if KB > 98
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 99
               MOVAPD 1456(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 99
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 99
               UNPCKHPD rb1, rb1
               MOVAPD 5952(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 99
               MOVAPD 5968(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 99
               MOVAPD 5984(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 99
               MOVAPD 6000(pA), rA3
            #endif
         #endif

         #if KB > 99
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 100
               MOVAPD 1472(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 100
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 100
               MOVAPD 6016(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 100
               MOVAPD 6032(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 100
               MOVAPD 6048(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 100
               MOVAPD 6064(pA), rA3
            #endif
         #endif
         #if KB > 100
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 101
               MOVAPD 1488(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 101
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 101
               UNPCKHPD rb1, rb1
               MOVAPD 6080(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 101
               MOVAPD 6096(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 101
               MOVAPD 6112(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 101
               MOVAPD 6128(pA), rA3
            #endif
         #endif

         #if KB > 101
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 102
               MOVAPD 1504(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 102
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 102
               MOVAPD 6144(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 102
               MOVAPD 6160(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 102
               MOVAPD 6176(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 102
               MOVAPD 6192(pA), rA3
            #endif
         #endif
         #if KB > 102
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 103
               MOVAPD 1520(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 103
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 103
               UNPCKHPD rb1, rb1
               MOVAPD 6208(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 103
               MOVAPD 6224(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 103
               MOVAPD 6240(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 103
               MOVAPD 6256(pA), rA3
            #endif
         #endif

         #if KB > 103
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 104
               MOVAPD 1536(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 104
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 104
               MOVAPD 6272(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 104
               MOVAPD 6288(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 104
               MOVAPD 6304(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 104
               MOVAPD 6320(pA), rA3
            #endif
         #endif
         #if KB > 104
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 105
               MOVAPD 1552(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 105
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 105
               UNPCKHPD rb1, rb1
               MOVAPD 6336(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 105
               MOVAPD 6352(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 105
               MOVAPD 6368(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 105
               MOVAPD 6384(pA), rA3
            #endif
         #endif

         #if KB > 105
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 106
               MOVAPD 1568(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 106
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 106
               MOVAPD 6400(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 106
               MOVAPD 6416(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 106
               MOVAPD 6432(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 106
               MOVAPD 6448(pA), rA3
            #endif
         #endif
         #if KB > 106
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 107
               MOVAPD 1584(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 107
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 107
               UNPCKHPD rb1, rb1
               MOVAPD 6464(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 107
               MOVAPD 6480(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 107
               MOVAPD 6496(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 107
               MOVAPD 6512(pA), rA3
            #endif
         #endif

         #if KB > 107
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 108
               MOVAPD 1600(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 108
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 108
               MOVAPD 6528(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 108
               MOVAPD 6544(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 108
               MOVAPD 6560(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 108
               MOVAPD 6576(pA), rA3
            #endif
         #endif
         #if KB > 108
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 109
               MOVAPD 1616(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 109
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 109
               UNPCKHPD rb1, rb1
               MOVAPD 6592(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 109
               MOVAPD 6608(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 109
               MOVAPD 6624(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 109
               MOVAPD 6640(pA), rA3
            #endif
         #endif

         #if KB > 109
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 110
               MOVAPD 1632(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 110
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 110
               MOVAPD 6656(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 110
               MOVAPD 6672(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 110
               MOVAPD 6688(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 110
               MOVAPD 6704(pA), rA3
            #endif
         #endif
         #if KB > 110
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 111
               MOVAPD 1648(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 111
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 111
               UNPCKHPD rb1, rb1
               MOVAPD 6720(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 111
               MOVAPD 6736(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 111
               MOVAPD 6752(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 111
               MOVAPD 6768(pA), rA3
            #endif
         #endif

         #if KB > 111
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 112
               MOVAPD 1664(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 112
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 112
               MOVAPD 6784(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 112
               MOVAPD 6800(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 112
               MOVAPD 6816(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 112
               MOVAPD 6832(pA), rA3
            #endif
         #endif
         #if KB > 112
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 113
               MOVAPD 1680(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 113
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 113
               UNPCKHPD rb1, rb1
               MOVAPD 6848(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 113
               MOVAPD 6864(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 113
               MOVAPD 6880(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 113
               MOVAPD 6896(pA), rA3
            #endif
         #endif

         #if KB > 113
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 114
               MOVAPD 1696(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 114
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 114
               MOVAPD 6912(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 114
               MOVAPD 6928(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 114
               MOVAPD 6944(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 114
               MOVAPD 6960(pA), rA3
            #endif
         #endif
         #if KB > 114
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 115
               MOVAPD 1712(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 115
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 115
               UNPCKHPD rb1, rb1
               MOVAPD 6976(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 115
               MOVAPD 6992(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 115
               MOVAPD 7008(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 115
               MOVAPD 7024(pA), rA3
            #endif
         #endif

         #if KB > 115
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 116
               MOVAPD 1728(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 116
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 116
               MOVAPD 7040(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 116
               MOVAPD 7056(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 116
               MOVAPD 7072(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 116
               MOVAPD 7088(pA), rA3
            #endif
         #endif
         #if KB > 116
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 117
               MOVAPD 1744(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 117
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 117
               UNPCKHPD rb1, rb1
               MOVAPD 7104(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 117
               MOVAPD 7120(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 117
               MOVAPD 7136(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 117
               MOVAPD 7152(pA), rA3
            #endif
         #endif

         #if KB > 117
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 118
               MOVAPD 1760(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 118
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 118
               MOVAPD 7168(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 118
               MOVAPD 7184(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 118
               MOVAPD 7200(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 118
               MOVAPD 7216(pA), rA3
            #endif
         #endif
         #if KB > 118
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 119
               MOVAPD 1776(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 119
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 119
               UNPCKHPD rb1, rb1
               MOVAPD 7232(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 119
               MOVAPD 7248(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 119
               MOVAPD 7264(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 119
               MOVAPD 7280(pA), rA3
            #endif
         #endif

         #if KB > 119
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 120
               MOVAPD 1792(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 120
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 120
               MOVAPD 7296(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 120
               MOVAPD 7312(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 120
               MOVAPD 7328(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 120
               MOVAPD 7344(pA), rA3
            #endif
         #endif
         #if KB > 120
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 121
               MOVAPD 1808(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 121
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 121
               UNPCKHPD rb1, rb1
               MOVAPD 7360(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 121
               MOVAPD 7376(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 121
               MOVAPD 7392(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 121
               MOVAPD 7408(pA), rA3
            #endif
         #endif

         #if KB > 121
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 122
               MOVAPD 1824(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 122
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 122
               MOVAPD 7424(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 122
               MOVAPD 7440(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 122
               MOVAPD 7456(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 122
               MOVAPD 7472(pA), rA3
            #endif
         #endif
         #if KB > 122
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 123
               MOVAPD 1840(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 123
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 123
               UNPCKHPD rb1, rb1
               MOVAPD 7488(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 123
               MOVAPD 7504(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 123
               MOVAPD 7520(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 123
               MOVAPD 7536(pA), rA3
            #endif
         #endif

         #if KB > 123
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 124
               MOVAPD 1856(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 124
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 124
               MOVAPD 7552(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 124
               MOVAPD 7568(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 124
               MOVAPD 7584(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 124
               MOVAPD 7600(pA), rA3
            #endif
         #endif
         #if KB > 124
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 125
               MOVAPD 1872(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 125
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 125
               UNPCKHPD rb1, rb1
               MOVAPD 7616(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 125
               MOVAPD 7632(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 125
               MOVAPD 7648(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 125
               MOVAPD 7664(pA), rA3
            #endif
         #endif

         #if KB > 125
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 126
               MOVAPD 1888(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 126
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 126
               MOVAPD 7680(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 126
               MOVAPD 7696(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 126
               MOVAPD 7712(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 126
               MOVAPD 7728(pA), rA3
            #endif
         #endif
         #if KB > 126
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 127
               MOVAPD 1904(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 127
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 127
               UNPCKHPD rb1, rb1
               MOVAPD 7744(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 127
               MOVAPD 7760(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 127
               MOVAPD 7776(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 127
               MOVAPD 7792(pA), rA3
            #endif
         #endif

         #if KB > 127
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 128
               MOVAPD 1920(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 128
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 128
               MOVAPD 7808(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 128
               MOVAPD 7824(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 128
               MOVAPD 7840(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 128
               MOVAPD 7856(pA), rA3
            #endif
         #endif
         #if KB > 128
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 129
               MOVAPD 1936(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 129
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 129
               UNPCKHPD rb1, rb1
               MOVAPD 7872(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 129
               MOVAPD 7888(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 129
               MOVAPD 7904(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 129
               MOVAPD 7920(pA), rA3
            #endif
         #endif

         #if KB > 129
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 130
               MOVAPD 1952(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 130
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 130
               MOVAPD 7936(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 130
               MOVAPD 7952(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 130
               MOVAPD 7968(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 130
               MOVAPD 7984(pA), rA3
            #endif
         #endif
         #if KB > 130
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 131
               MOVAPD 1968(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 131
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 131
               UNPCKHPD rb1, rb1
               MOVAPD 8000(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 131
               MOVAPD 8016(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 131
               MOVAPD 8032(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 131
               MOVAPD 8048(pA), rA3
            #endif
         #endif

         #if KB > 131
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 132
               MOVAPD 1984(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 132
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 132
               MOVAPD 8064(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 132
               MOVAPD 8080(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 132
               MOVAPD 8096(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 132
               MOVAPD 8112(pA), rA3
            #endif
         #endif
         #if KB > 132
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 133
               MOVAPD 2000(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 133
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 133
               UNPCKHPD rb1, rb1
               MOVAPD 8128(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 133
               MOVAPD 8144(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 133
               MOVAPD 8160(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 133
               MOVAPD 8176(pA), rA3
            #endif
         #endif

         #if KB > 133
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 134
               MOVAPD 2016(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 134
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 134
               MOVAPD 8192(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 134
               MOVAPD 8208(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 134
               MOVAPD 8224(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 134
               MOVAPD 8240(pA), rA3
            #endif
         #endif
         #if KB > 134
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 135
               MOVAPD 2032(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 135
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 135
               UNPCKHPD rb1, rb1
               MOVAPD 8256(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 135
               MOVAPD 8272(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 135
               MOVAPD 8288(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 135
               MOVAPD 8304(pA), rA3
            #endif
         #endif

         #if KB > 135
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 136
               MOVAPD 2048(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 136
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 136
               MOVAPD 8320(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 136
               MOVAPD 8336(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 136
               MOVAPD 8352(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 136
               MOVAPD 8368(pA), rA3
            #endif
         #endif
         #if KB > 136
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 137
               MOVAPD 2064(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 137
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 137
               UNPCKHPD rb1, rb1
               MOVAPD 8384(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 137
               MOVAPD 8400(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 137
               MOVAPD 8416(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 137
               MOVAPD 8432(pA), rA3
            #endif
         #endif

         #if KB > 137
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 138
               MOVAPD 2080(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 138
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 138
               MOVAPD 8448(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 138
               MOVAPD 8464(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 138
               MOVAPD 8480(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 138
               MOVAPD 8496(pA), rA3
            #endif
         #endif
         #if KB > 138
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 139
               MOVAPD 2096(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 139
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 139
               UNPCKHPD rb1, rb1
               MOVAPD 8512(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 139
               MOVAPD 8528(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 139
               MOVAPD 8544(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 139
               MOVAPD 8560(pA), rA3
            #endif
         #endif

         #if KB > 139
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 140
               MOVAPD 2112(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 140
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 140
               MOVAPD 8576(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 140
               MOVAPD 8592(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 140
               MOVAPD 8608(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 140
               MOVAPD 8624(pA), rA3
            #endif
         #endif
         #if KB > 140
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 141
               MOVAPD 2128(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 141
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 141
               UNPCKHPD rb1, rb1
               MOVAPD 8640(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 141
               MOVAPD 8656(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 141
               MOVAPD 8672(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 141
               MOVAPD 8688(pA), rA3
            #endif
         #endif

         #if KB > 141
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 142
               MOVAPD 2144(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 142
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 142
               MOVAPD 8704(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 142
               MOVAPD 8720(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 142
               MOVAPD 8736(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 142
               MOVAPD 8752(pA), rA3
            #endif
         #endif
         #if KB > 142
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 143
               MOVAPD 2160(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 143
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 143
               UNPCKHPD rb1, rb1
               MOVAPD 8768(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 143
               MOVAPD 8784(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 143
               MOVAPD 8800(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 143
               MOVAPD 8816(pA), rA3
            #endif
         #endif

         #if KB > 143
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 144
               MOVAPD 2176(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 144
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 144
               MOVAPD 8832(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 144
               MOVAPD 8848(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 144
               MOVAPD 8864(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 144
               MOVAPD 8880(pA), rA3
            #endif
         #endif
         #if KB > 144
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 145
               MOVAPD 2192(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 145
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 145
               UNPCKHPD rb1, rb1
               MOVAPD 8896(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 145
               MOVAPD 8912(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 145
               MOVAPD 8928(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 145
               MOVAPD 8944(pA), rA3
            #endif
         #endif

         #if KB > 145
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 146
               MOVAPD 2208(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 146
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 146
               MOVAPD 8960(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 146
               MOVAPD 8976(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 146
               MOVAPD 8992(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 146
               MOVAPD 9008(pA), rA3
            #endif
         #endif
         #if KB > 146
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 147
               MOVAPD 2224(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 147
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 147
               UNPCKHPD rb1, rb1
               MOVAPD 9024(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 147
               MOVAPD 9040(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 147
               MOVAPD 9056(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 147
               MOVAPD 9072(pA), rA3
            #endif
         #endif

         #if KB > 147
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 148
               MOVAPD 2240(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 148
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 148
               MOVAPD 9088(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 148
               MOVAPD 9104(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 148
               MOVAPD 9120(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 148
               MOVAPD 9136(pA), rA3
            #endif
         #endif
         #if KB > 148
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 149
               MOVAPD 2256(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 149
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 149
               UNPCKHPD rb1, rb1
               MOVAPD 9152(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 149
               MOVAPD 9168(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 149
               MOVAPD 9184(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 149
               MOVAPD 9200(pA), rA3
            #endif
         #endif

         #if KB > 149
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 150
               MOVAPD 2272(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 150
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 150
               MOVAPD 9216(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 150
               MOVAPD 9232(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 150
               MOVAPD 9248(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 150
               MOVAPD 9264(pA), rA3
            #endif
         #endif
         #if KB > 150
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 151
               MOVAPD 2288(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 151
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 151
               UNPCKHPD rb1, rb1
               MOVAPD 9280(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 151
               MOVAPD 9296(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 151
               MOVAPD 9312(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 151
               MOVAPD 9328(pA), rA3
            #endif
         #endif

         #if KB > 151
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 152
               MOVAPD 2304(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 152
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 152
               MOVAPD 9344(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 152
               MOVAPD 9360(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 152
               MOVAPD 9376(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 152
               MOVAPD 9392(pA), rA3
            #endif
         #endif
         #if KB > 152
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 153
               MOVAPD 2320(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 153
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 153
               UNPCKHPD rb1, rb1
               MOVAPD 9408(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 153
               MOVAPD 9424(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 153
               MOVAPD 9440(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 153
               MOVAPD 9456(pA), rA3
            #endif
         #endif

         #if KB > 153
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 154
               MOVAPD 2336(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 154
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 154
               MOVAPD 9472(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 154
               MOVAPD 9488(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 154
               MOVAPD 9504(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 154
               MOVAPD 9520(pA), rA3
            #endif
         #endif
         #if KB > 154
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 155
               MOVAPD 2352(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 155
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 155
               UNPCKHPD rb1, rb1
               MOVAPD 9536(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 155
               MOVAPD 9552(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 155
               MOVAPD 9568(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 155
               MOVAPD 9584(pA), rA3
            #endif
         #endif

         #if KB > 155
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 156
               MOVAPD 2368(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 156
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 156
               MOVAPD 9600(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 156
               MOVAPD 9616(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 156
               MOVAPD 9632(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 156
               MOVAPD 9648(pA), rA3
            #endif
         #endif
         #if KB > 156
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 157
               MOVAPD 2384(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 157
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 157
               UNPCKHPD rb1, rb1
               MOVAPD 9664(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 157
               MOVAPD 9680(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 157
               MOVAPD 9696(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 157
               MOVAPD 9712(pA), rA3
            #endif
         #endif

         #if KB > 157
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 158
               MOVAPD 2400(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 158
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 158
               MOVAPD 9728(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 158
               MOVAPD 9744(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 158
               MOVAPD 9760(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 158
               MOVAPD 9776(pA), rA3
            #endif
         #endif
         #if KB > 158
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 159
               MOVAPD 2416(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 159
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 159
               UNPCKHPD rb1, rb1
               MOVAPD 9792(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 159
               MOVAPD 9808(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 159
               MOVAPD 9824(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 159
               MOVAPD 9840(pA), rA3
            #endif
         #endif

         #if KB > 159
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 160
               MOVAPD 2432(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 160
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 160
               MOVAPD 9856(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 160
               MOVAPD 9872(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 160
               MOVAPD 9888(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 160
               MOVAPD 9904(pA), rA3
            #endif
         #endif
         #if KB > 160
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 161
               MOVAPD 2448(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 161
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 161
               UNPCKHPD rb1, rb1
               MOVAPD 9920(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 161
               MOVAPD 9936(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 161
               MOVAPD 9952(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 161
               MOVAPD 9968(pA), rA3
            #endif
         #endif

         #if KB > 161
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 162
               MOVAPD 2464(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 162
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 162
               MOVAPD 9984(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 162
               MOVAPD 10000(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 162
               MOVAPD 10016(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 162
               MOVAPD 10032(pA), rA3
            #endif
         #endif
         #if KB > 162
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 163
               MOVAPD 2480(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 163
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 163
               UNPCKHPD rb1, rb1
               MOVAPD 10048(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 163
               MOVAPD 10064(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 163
               MOVAPD 10080(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 163
               MOVAPD 10096(pA), rA3
            #endif
         #endif

         #if KB > 163
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 164
               MOVAPD 2496(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 164
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 164
               MOVAPD 10112(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 164
               MOVAPD 10128(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 164
               MOVAPD 10144(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 164
               MOVAPD 10160(pA), rA3
            #endif
         #endif
         #if KB > 164
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 165
               MOVAPD 2512(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 165
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 165
               UNPCKHPD rb1, rb1
               MOVAPD 10176(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 165
               MOVAPD 10192(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 165
               MOVAPD 10208(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 165
               MOVAPD 10224(pA), rA3
            #endif
         #endif

         #if KB > 165
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 166
               MOVAPD 2528(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 166
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 166
               MOVAPD 10240(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 166
               MOVAPD 10256(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 166
               MOVAPD 10272(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 166
               MOVAPD 10288(pA), rA3
            #endif
         #endif
         #if KB > 166
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 167
               MOVAPD 2544(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 167
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 167
               UNPCKHPD rb1, rb1
               MOVAPD 10304(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 167
               MOVAPD 10320(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 167
               MOVAPD 10336(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 167
               MOVAPD 10352(pA), rA3
            #endif
         #endif

         #if KB > 167
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 168
               MOVAPD 2560(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 168
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 168
               MOVAPD 10368(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 168
               MOVAPD 10384(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 168
               MOVAPD 10400(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 168
               MOVAPD 10416(pA), rA3
            #endif
         #endif
         #if KB > 168
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 169
               MOVAPD 2576(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 169
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 169
               UNPCKHPD rb1, rb1
               MOVAPD 10432(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 169
               MOVAPD 10448(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 169
               MOVAPD 10464(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 169
               MOVAPD 10480(pA), rA3
            #endif
         #endif

         #if KB > 169
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 170
               MOVAPD 2592(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 170
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 170
               MOVAPD 10496(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 170
               MOVAPD 10512(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 170
               MOVAPD 10528(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 170
               MOVAPD 10544(pA), rA3
            #endif
         #endif
         #if KB > 170
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 171
               MOVAPD 2608(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 171
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 171
               UNPCKHPD rb1, rb1
               MOVAPD 10560(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 171
               MOVAPD 10576(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 171
               MOVAPD 10592(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 171
               MOVAPD 10608(pA), rA3
            #endif
         #endif

         #if KB > 171
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 172
               MOVAPD 2624(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 172
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 172
               MOVAPD 10624(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 172
               MOVAPD 10640(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 172
               MOVAPD 10656(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 172
               MOVAPD 10672(pA), rA3
            #endif
         #endif
         #if KB > 172
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 173
               MOVAPD 2640(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 173
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 173
               UNPCKHPD rb1, rb1
               MOVAPD 10688(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 173
               MOVAPD 10704(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 173
               MOVAPD 10720(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 173
               MOVAPD 10736(pA), rA3
            #endif
         #endif

         #if KB > 173
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 174
               MOVAPD 2656(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 174
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 174
               MOVAPD 10752(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 174
               MOVAPD 10768(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 174
               MOVAPD 10784(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 174
               MOVAPD 10800(pA), rA3
            #endif
         #endif
         #if KB > 174
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 175
               MOVAPD 2672(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 175
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 175
               UNPCKHPD rb1, rb1
               MOVAPD 10816(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 175
               MOVAPD 10832(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 175
               MOVAPD 10848(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 175
               MOVAPD 10864(pA), rA3
            #endif
         #endif

         #if KB > 175
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 176
               MOVAPD 2688(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 176
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 176
               MOVAPD 10880(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 176
               MOVAPD 10896(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 176
               MOVAPD 10912(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 176
               MOVAPD 10928(pA), rA3
            #endif
         #endif
         #if KB > 176
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 177
               MOVAPD 2704(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 177
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 177
               UNPCKHPD rb1, rb1
               MOVAPD 10944(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 177
               MOVAPD 10960(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 177
               MOVAPD 10976(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 177
               MOVAPD 10992(pA), rA3
            #endif
         #endif

         #if KB > 177
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 178
               MOVAPD 2720(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 178
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 178
               MOVAPD 11008(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 178
               MOVAPD 11024(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 178
               MOVAPD 11040(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 178
               MOVAPD 11056(pA), rA3
            #endif
         #endif
         #if KB > 178
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 179
               MOVAPD 2736(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 179
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 179
               UNPCKHPD rb1, rb1
               MOVAPD 11072(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 179
               MOVAPD 11088(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 179
               MOVAPD 11104(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 179
               MOVAPD 11120(pA), rA3
            #endif
         #endif

         #if KB > 179
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 180
               MOVAPD 2752(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 180
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 180
               MOVAPD 11136(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 180
               MOVAPD 11152(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 180
               MOVAPD 11168(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 180
               MOVAPD 11184(pA), rA3
            #endif
         #endif
         #if KB > 180
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 181
               MOVAPD 2768(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 181
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 181
               UNPCKHPD rb1, rb1
               MOVAPD 11200(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 181
               MOVAPD 11216(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 181
               MOVAPD 11232(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 181
               MOVAPD 11248(pA), rA3
            #endif
         #endif

         #if KB > 181
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 182
               MOVAPD 2784(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 182
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 182
               MOVAPD 11264(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 182
               MOVAPD 11280(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 182
               MOVAPD 11296(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 182
               MOVAPD 11312(pA), rA3
            #endif
         #endif
         #if KB > 182
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 183
               MOVAPD 2800(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 183
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 183
               UNPCKHPD rb1, rb1
               MOVAPD 11328(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 183
               MOVAPD 11344(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 183
               MOVAPD 11360(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 183
               MOVAPD 11376(pA), rA3
            #endif
         #endif

         #if KB > 183
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 184
               MOVAPD 2816(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 184
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 184
               MOVAPD 11392(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 184
               MOVAPD 11408(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 184
               MOVAPD 11424(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 184
               MOVAPD 11440(pA), rA3
            #endif
         #endif
         #if KB > 184
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 185
               MOVAPD 2832(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 185
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 185
               UNPCKHPD rb1, rb1
               MOVAPD 11456(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 185
               MOVAPD 11472(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 185
               MOVAPD 11488(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 185
               MOVAPD 11504(pA), rA3
            #endif
         #endif

         #if KB > 185
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 186
               MOVAPD 2848(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 186
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 186
               MOVAPD 11520(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 186
               MOVAPD 11536(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 186
               MOVAPD 11552(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 186
               MOVAPD 11568(pA), rA3
            #endif
         #endif
         #if KB > 186
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 187
               MOVAPD 2864(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 187
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 187
               UNPCKHPD rb1, rb1
               MOVAPD 11584(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 187
               MOVAPD 11600(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 187
               MOVAPD 11616(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 187
               MOVAPD 11632(pA), rA3
            #endif
         #endif

         #if KB > 187
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 188
               MOVAPD 2880(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 188
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 188
               MOVAPD 11648(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 188
               MOVAPD 11664(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 188
               MOVAPD 11680(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 188
               MOVAPD 11696(pA), rA3
            #endif
         #endif
         #if KB > 188
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 189
               MOVAPD 2896(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 189
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 189
               UNPCKHPD rb1, rb1
               MOVAPD 11712(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 189
               MOVAPD 11728(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 189
               MOVAPD 11744(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 189
               MOVAPD 11760(pA), rA3
            #endif
         #endif

         #if KB > 189
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 190
               MOVAPD 2912(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 190
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 190
               MOVAPD 11776(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 190
               MOVAPD 11792(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 190
               MOVAPD 11808(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 190
               MOVAPD 11824(pA), rA3
            #endif
         #endif
         #if KB > 190
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 191
               MOVAPD 2928(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 191
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 191
               UNPCKHPD rb1, rb1
               MOVAPD 11840(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 191
               MOVAPD 11856(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 191
               MOVAPD 11872(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 191
               MOVAPD 11888(pA), rA3
            #endif
         #endif

         #if KB > 191
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 192
               MOVAPD 2944(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 192
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 192
               MOVAPD 11904(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 192
               MOVAPD 11920(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 192
               MOVAPD 11936(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 192
               MOVAPD 11952(pA), rA3
            #endif
         #endif
         #if KB > 192
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 193
               MOVAPD 2960(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 193
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 193
               UNPCKHPD rb1, rb1
               MOVAPD 11968(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 193
               MOVAPD 11984(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 193
               MOVAPD 12000(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 193
               MOVAPD 12016(pA), rA3
            #endif
         #endif

         #if KB > 193
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 194
               MOVAPD 2976(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 194
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 194
               MOVAPD 12032(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 194
               MOVAPD 12048(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 194
               MOVAPD 12064(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 194
               MOVAPD 12080(pA), rA3
            #endif
         #endif
         #if KB > 194
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 195
               MOVAPD 2992(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 195
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 195
               UNPCKHPD rb1, rb1
               MOVAPD 12096(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 195
               MOVAPD 12112(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 195
               MOVAPD 12128(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 195
               MOVAPD 12144(pA), rA3
            #endif
         #endif

         #if KB > 195
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 196
               MOVAPD 3008(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 196
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 196
               MOVAPD 12160(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 196
               MOVAPD 12176(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 196
               MOVAPD 12192(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 196
               MOVAPD 12208(pA), rA3
            #endif
         #endif
         #if KB > 196
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 197
               MOVAPD 3024(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 197
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 197
               UNPCKHPD rb1, rb1
               MOVAPD 12224(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 197
               MOVAPD 12240(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 197
               MOVAPD 12256(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 197
               MOVAPD 12272(pA), rA3
            #endif
         #endif

         #if KB > 197
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 198
               MOVAPD 3040(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 198
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 198
               MOVAPD 12288(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 198
               MOVAPD 12304(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 198
               MOVAPD 12320(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 198
               MOVAPD 12336(pA), rA3
            #endif
         #endif
         #if KB > 198
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 199
               MOVAPD 3056(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 199
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 199
               UNPCKHPD rb1, rb1
               MOVAPD 12352(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 199
               MOVAPD 12368(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 199
               MOVAPD 12384(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 199
               MOVAPD 12400(pA), rA3
            #endif
         #endif

         #if KB > 199
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 200
               MOVAPD 3072(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 200
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 200
               MOVAPD 12416(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 200
               MOVAPD 12432(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 200
               MOVAPD 12448(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 200
               MOVAPD 12464(pA), rA3
            #endif
         #endif
         #if KB > 200
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 201
               MOVAPD 3088(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 201
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 201
               UNPCKHPD rb1, rb1
               MOVAPD 12480(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 201
               MOVAPD 12496(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 201
               MOVAPD 12512(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 201
               MOVAPD 12528(pA), rA3
            #endif
         #endif

         #if KB > 201
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 202
               MOVAPD 3104(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 202
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 202
               MOVAPD 12544(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 202
               MOVAPD 12560(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 202
               MOVAPD 12576(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 202
               MOVAPD 12592(pA), rA3
            #endif
         #endif
         #if KB > 202
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 203
               MOVAPD 3120(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 203
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 203
               UNPCKHPD rb1, rb1
               MOVAPD 12608(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 203
               MOVAPD 12624(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 203
               MOVAPD 12640(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 203
               MOVAPD 12656(pA), rA3
            #endif
         #endif

         #if KB > 203
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 204
               MOVAPD 3136(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 204
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 204
               MOVAPD 12672(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 204
               MOVAPD 12688(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 204
               MOVAPD 12704(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 204
               MOVAPD 12720(pA), rA3
            #endif
         #endif
         #if KB > 204
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 205
               MOVAPD 3152(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 205
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 205
               UNPCKHPD rb1, rb1
               MOVAPD 12736(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 205
               MOVAPD 12752(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 205
               MOVAPD 12768(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 205
               MOVAPD 12784(pA), rA3
            #endif
         #endif

         #if KB > 205
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 206
               MOVAPD 3168(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 206
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 206
               MOVAPD 12800(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 206
               MOVAPD 12816(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 206
               MOVAPD 12832(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 206
               MOVAPD 12848(pA), rA3
            #endif
         #endif
         #if KB > 206
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 207
               MOVAPD 3184(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 207
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 207
               UNPCKHPD rb1, rb1
               MOVAPD 12864(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 207
               MOVAPD 12880(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 207
               MOVAPD 12896(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 207
               MOVAPD 12912(pA), rA3
            #endif
         #endif

         #if KB > 207
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 208
               MOVAPD 3200(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 208
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 208
               MOVAPD 12928(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 208
               MOVAPD 12944(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 208
               MOVAPD 12960(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 208
               MOVAPD 12976(pA), rA3
            #endif
         #endif
         #if KB > 208
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 209
               MOVAPD 3216(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 209
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 209
               UNPCKHPD rb1, rb1
               MOVAPD 12992(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 209
               MOVAPD 13008(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 209
               MOVAPD 13024(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 209
               MOVAPD 13040(pA), rA3
            #endif
         #endif

         #if KB > 209
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 210
               MOVAPD 3232(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 210
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 210
               MOVAPD 13056(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 210
               MOVAPD 13072(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 210
               MOVAPD 13088(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 210
               MOVAPD 13104(pA), rA3
            #endif
         #endif
         #if KB > 210
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 211
               MOVAPD 3248(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 211
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 211
               UNPCKHPD rb1, rb1
               MOVAPD 13120(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 211
               MOVAPD 13136(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 211
               MOVAPD 13152(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 211
               MOVAPD 13168(pA), rA3
            #endif
         #endif

         #if KB > 211
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 212
               MOVAPD 3264(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 212
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 212
               MOVAPD 13184(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 212
               MOVAPD 13200(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 212
               MOVAPD 13216(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 212
               MOVAPD 13232(pA), rA3
            #endif
         #endif
         #if KB > 212
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 213
               MOVAPD 3280(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 213
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 213
               UNPCKHPD rb1, rb1
               MOVAPD 13248(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 213
               MOVAPD 13264(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 213
               MOVAPD 13280(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 213
               MOVAPD 13296(pA), rA3
            #endif
         #endif

         #if KB > 213
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 214
               MOVAPD 3296(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 214
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 214
               MOVAPD 13312(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 214
               MOVAPD 13328(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 214
               MOVAPD 13344(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 214
               MOVAPD 13360(pA), rA3
            #endif
         #endif
         #if KB > 214
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 215
               MOVAPD 3312(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 215
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 215
               UNPCKHPD rb1, rb1
               MOVAPD 13376(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 215
               MOVAPD 13392(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 215
               MOVAPD 13408(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 215
               MOVAPD 13424(pA), rA3
            #endif
         #endif

         #if KB > 215
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 216
               MOVAPD 3328(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 216
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 216
               MOVAPD 13440(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 216
               MOVAPD 13456(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 216
               MOVAPD 13472(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 216
               MOVAPD 13488(pA), rA3
            #endif
         #endif
         #if KB > 216
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 217
               MOVAPD 3344(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 217
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 217
               UNPCKHPD rb1, rb1
               MOVAPD 13504(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 217
               MOVAPD 13520(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 217
               MOVAPD 13536(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 217
               MOVAPD 13552(pA), rA3
            #endif
         #endif

         #if KB > 217
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 218
               MOVAPD 3360(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 218
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 218
               MOVAPD 13568(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 218
               MOVAPD 13584(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 218
               MOVAPD 13600(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 218
               MOVAPD 13616(pA), rA3
            #endif
         #endif
         #if KB > 218
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 219
               MOVAPD 3376(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 219
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 219
               UNPCKHPD rb1, rb1
               MOVAPD 13632(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 219
               MOVAPD 13648(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 219
               MOVAPD 13664(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 219
               MOVAPD 13680(pA), rA3
            #endif
         #endif

         #if KB > 219
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 220
               MOVAPD 3392(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 220
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 220
               MOVAPD 13696(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 220
               MOVAPD 13712(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 220
               MOVAPD 13728(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 220
               MOVAPD 13744(pA), rA3
            #endif
         #endif
         #if KB > 220
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 221
               MOVAPD 3408(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 221
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 221
               UNPCKHPD rb1, rb1
               MOVAPD 13760(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 221
               MOVAPD 13776(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 221
               MOVAPD 13792(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 221
               MOVAPD 13808(pA), rA3
            #endif
         #endif

         #if KB > 221
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 222
               MOVAPD 3424(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 222
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 222
               MOVAPD 13824(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 222
               MOVAPD 13840(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 222
               MOVAPD 13856(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 222
               MOVAPD 13872(pA), rA3
            #endif
         #endif
         #if KB > 222
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 223
               MOVAPD 3440(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 223
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 223
               UNPCKHPD rb1, rb1
               MOVAPD 13888(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 223
               MOVAPD 13904(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 223
               MOVAPD 13920(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 223
               MOVAPD 13936(pA), rA3
            #endif
         #endif

         #if KB > 223
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 224
               MOVAPD 3456(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 224
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 224
               MOVAPD 13952(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 224
               MOVAPD 13968(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 224
               MOVAPD 13984(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 224
               MOVAPD 14000(pA), rA3
            #endif
         #endif
         #if KB > 224
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 225
               MOVAPD 3472(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 225
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 225
               UNPCKHPD rb1, rb1
               MOVAPD 14016(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 225
               MOVAPD 14032(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 225
               MOVAPD 14048(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 225
               MOVAPD 14064(pA), rA3
            #endif
         #endif

         #if KB > 225
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 226
               MOVAPD 3488(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 226
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 226
               MOVAPD 14080(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 226
               MOVAPD 14096(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 226
               MOVAPD 14112(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 226
               MOVAPD 14128(pA), rA3
            #endif
         #endif
         #if KB > 226
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 227
               MOVAPD 3504(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 227
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 227
               UNPCKHPD rb1, rb1
               MOVAPD 14144(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 227
               MOVAPD 14160(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 227
               MOVAPD 14176(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 227
               MOVAPD 14192(pA), rA3
            #endif
         #endif

         #if KB > 227
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 228
               MOVAPD 3520(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 228
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 228
               MOVAPD 14208(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 228
               MOVAPD 14224(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 228
               MOVAPD 14240(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 228
               MOVAPD 14256(pA), rA3
            #endif
         #endif
         #if KB > 228
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 229
               MOVAPD 3536(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 229
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 229
               UNPCKHPD rb1, rb1
               MOVAPD 14272(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 229
               MOVAPD 14288(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 229
               MOVAPD 14304(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 229
               MOVAPD 14320(pA), rA3
            #endif
         #endif

         #if KB > 229
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 230
               MOVAPD 3552(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 230
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 230
               MOVAPD 14336(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 230
               MOVAPD 14352(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 230
               MOVAPD 14368(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 230
               MOVAPD 14384(pA), rA3
            #endif
         #endif
         #if KB > 230
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 231
               MOVAPD 3568(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 231
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 231
               UNPCKHPD rb1, rb1
               MOVAPD 14400(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 231
               MOVAPD 14416(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 231
               MOVAPD 14432(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 231
               MOVAPD 14448(pA), rA3
            #endif
         #endif

         #if KB > 231
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 232
               MOVAPD 3584(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 232
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 232
               MOVAPD 14464(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 232
               MOVAPD 14480(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 232
               MOVAPD 14496(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 232
               MOVAPD 14512(pA), rA3
            #endif
         #endif
         #if KB > 232
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 233
               MOVAPD 3600(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 233
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 233
               UNPCKHPD rb1, rb1
               MOVAPD 14528(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 233
               MOVAPD 14544(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 233
               MOVAPD 14560(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 233
               MOVAPD 14576(pA), rA3
            #endif
         #endif

         #if KB > 233
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 234
               MOVAPD 3616(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 234
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 234
               MOVAPD 14592(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 234
               MOVAPD 14608(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 234
               MOVAPD 14624(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 234
               MOVAPD 14640(pA), rA3
            #endif
         #endif
         #if KB > 234
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 235
               MOVAPD 3632(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 235
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 235
               UNPCKHPD rb1, rb1
               MOVAPD 14656(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 235
               MOVAPD 14672(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 235
               MOVAPD 14688(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 235
               MOVAPD 14704(pA), rA3
            #endif
         #endif

         #if KB > 235
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 236
               MOVAPD 3648(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 236
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 236
               MOVAPD 14720(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 236
               MOVAPD 14736(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 236
               MOVAPD 14752(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 236
               MOVAPD 14768(pA), rA3
            #endif
         #endif
         #if KB > 236
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 237
               MOVAPD 3664(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 237
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 237
               UNPCKHPD rb1, rb1
               MOVAPD 14784(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 237
               MOVAPD 14800(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 237
               MOVAPD 14816(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 237
               MOVAPD 14832(pA), rA3
            #endif
         #endif

         #if KB > 237
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 238
               MOVAPD 3680(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 238
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 238
               MOVAPD 14848(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 238
               MOVAPD 14864(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 238
               MOVAPD 14880(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 238
               MOVAPD 14896(pA), rA3
            #endif
         #endif
         #if KB > 238
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 239
               MOVAPD 3696(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 239
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 239
               UNPCKHPD rb1, rb1
               MOVAPD 14912(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 239
               MOVAPD 14928(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 239
               MOVAPD 14944(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 239
               MOVAPD 14960(pA), rA3
            #endif
         #endif

         #if KB > 239
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 240
               MOVAPD 3712(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 240
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 240
               MOVAPD 14976(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 240
               MOVAPD 14992(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 240
               MOVAPD 15008(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 240
               MOVAPD 15024(pA), rA3
            #endif
         #endif
         #if KB > 240
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 241
               MOVAPD 3728(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 241
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 241
               UNPCKHPD rb1, rb1
               MOVAPD 15040(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 241
               MOVAPD 15056(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 241
               MOVAPD 15072(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 241
               MOVAPD 15088(pA), rA3
            #endif
         #endif

         #if KB > 241
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 242
               MOVAPD 3744(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 242
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 242
               MOVAPD 15104(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 242
               MOVAPD 15120(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 242
               MOVAPD 15136(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 242
               MOVAPD 15152(pA), rA3
            #endif
         #endif
         #if KB > 242
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 243
               MOVAPD 3760(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 243
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 243
               UNPCKHPD rb1, rb1
               MOVAPD 15168(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 243
               MOVAPD 15184(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 243
               MOVAPD 15200(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 243
               MOVAPD 15216(pA), rA3
            #endif
         #endif

         #if KB > 243
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 244
               MOVAPD 3776(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 244
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 244
               MOVAPD 15232(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 244
               MOVAPD 15248(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 244
               MOVAPD 15264(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 244
               MOVAPD 15280(pA), rA3
            #endif
         #endif
         #if KB > 244
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 245
               MOVAPD 3792(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 245
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 245
               UNPCKHPD rb1, rb1
               MOVAPD 15296(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 245
               MOVAPD 15312(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 245
               MOVAPD 15328(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 245
               MOVAPD 15344(pA), rA3
            #endif
         #endif

         #if KB > 245
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 246
               MOVAPD 3808(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 246
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 246
               MOVAPD 15360(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 246
               MOVAPD 15376(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 246
               MOVAPD 15392(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 246
               MOVAPD 15408(pA), rA3
            #endif
         #endif
         #if KB > 246
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 247
               MOVAPD 3824(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 247
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 247
               UNPCKHPD rb1, rb1
               MOVAPD 15424(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 247
               MOVAPD 15440(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 247
               MOVAPD 15456(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 247
               MOVAPD 15472(pA), rA3
            #endif
         #endif

         #if KB > 247
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 248
               MOVAPD 3840(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 248
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 248
               MOVAPD 15488(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 248
               MOVAPD 15504(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 248
               MOVAPD 15520(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 248
               MOVAPD 15536(pA), rA3
            #endif
         #endif
         #if KB > 248
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 249
               MOVAPD 3856(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 249
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 249
               UNPCKHPD rb1, rb1
               MOVAPD 15552(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 249
               MOVAPD 15568(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 249
               MOVAPD 15584(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 249
               MOVAPD 15600(pA), rA3
            #endif
         #endif

         #if KB > 249
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 250
               MOVAPD 3872(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 250
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 250
               MOVAPD 15616(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 250
               MOVAPD 15632(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 250
               MOVAPD 15648(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 250
               MOVAPD 15664(pA), rA3
            #endif
         #endif
         #if KB > 250
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 251
               MOVAPD 3888(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 251
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 251
               UNPCKHPD rb1, rb1
               MOVAPD 15680(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 251
               MOVAPD 15696(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 251
               MOVAPD 15712(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 251
               MOVAPD 15728(pA), rA3
            #endif
         #endif

         #if KB > 251
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 252
               MOVAPD 3904(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 252
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 252
               MOVAPD 15744(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 252
               MOVAPD 15760(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 252
               MOVAPD 15776(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 252
               MOVAPD 15792(pA), rA3
            #endif
         #endif
         #if KB > 252
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 253
               MOVAPD 3920(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 253
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 253
               UNPCKHPD rb1, rb1
               MOVAPD 15808(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 253
               MOVAPD 15824(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 253
               MOVAPD 15840(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 253
               MOVAPD 15856(pA), rA3
            #endif
         #endif

         #if KB > 253
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 254
               MOVAPD 3936(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 254
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 254
               MOVAPD 15872(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 254
               MOVAPD 15888(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 254
               MOVAPD 15904(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 254
               MOVAPD 15920(pA), rA3
            #endif
         #endif
         #if KB > 254
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 255
               MOVAPD 3952(pB), rb1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 255
               movddup rb1, rB0
            #endif

            mulpd rB1, rA0
            addpd rA0, rC01
            #if KB > 255
               UNPCKHPD rb1, rb1
               MOVAPD 15936(pA), rA0
            #endif
            mulpd rB1, rA1
            addpd rA1, rC11
            #if KB > 255
               MOVAPD 15952(pA), rA1
            #endif
            mulpd rB1, rA2
            addpd rA2, rC21
            #if KB > 255
               MOVAPD 15968(pA), rA2
            #endif
            mulpd rA3, rB1
            addpd rB1, rC31
            #if KB > 255
               MOVAPD 15984(pA), rA3
            #endif
         #endif

         #if KB > 255
            MOVAPD rA0, m0
            mulpd rB0, m0
            addpd m0, rC00
            #if KB > 256
               MOVAPD 3968(pB), rB1
            #endif
            MOVAPD rA1, m0
            mulpd rB0, m0
            addpd m0, rC10
            MOVAPD rA2, m0
            mulpd rB0, m0
            addpd m0, rC20
            UNPCKHPD rb1, rb1
            mulpd rA3, rB0
            addpd rB0, rC30
            #if KB > 256
               movddup rB1, rB0
            #endif

            mulpd rb1, rA0
            addpd rA0, rC01
            #if KB > 256
               MOVAPD 16000(pA), rA0
               UNPCKHPD rB1, rB1
            #endif
            mulpd rb1, rA1
            addpd rA1, rC11
            #if KB > 256
               MOVAPD 16016(pA), rA1
            #endif
            mulpd rb1, rA2
            addpd rA2, rC21
            #if KB > 256
               MOVAPD 16032(pA), rA2
            #endif
            mulpd rA3, rb1
            addpd rb1, rC31
            #if KB > 256
               MOVAPD 16048(pA), rA3
            #endif
         #endif
/*
 *       Write answer out to C
 */
         #ifdef BETCOP
            BETCOP (pC), rC00
         #endif
         MOVAPD rC00, (pC)
         #ifdef BETCOP
            BETCOP 16(pC), rC10
         #endif
         MOVAPD rC10, 16(pC)
         #ifdef BETCOP
            BETCOP 32(pC), rC20
         #endif
         MOVAPD rC20, 32(pC)
         #ifdef BETCOP
            BETCOP 48(pC), rC30
         #endif
         MOVAPD rC30, 48(pC)
         #ifdef BETCOP
            BETCOP 64(pC), rC01
         #endif
         MOVAPD rC01, 64(pC)
         #ifdef BETCOP
            BETCOP 80(pC), rC11
         #endif
         MOVAPD rC11, 80(pC)
         #ifdef BETCOP
            BETCOP 96(pC), rC21
         #endif
         MOVAPD rC21, 96(pC)
         #ifdef BETCOP
            BETCOP 112(pC), rC31
         #endif
         MOVAPD rC31, 112(pC)
         add incPF, pC
         mov pA0, pA
         add incB, pB
      dec nnu
      jnz MNLOOP

      mov nnu0, nnu
      add incA, pA0
      mov pA0, pA
      mov pB0, pB
   dec nmu
   jnz MNLOOP

/* DONE: */
   movq    (%rsp), %rbp
   movq    8(%rsp), %rbx
   movq    16(%rsp), %r12
   movq    24(%rsp), %r13
   movq    32(%rsp), %r14
   add $FSIZE, %rsp
   ret
