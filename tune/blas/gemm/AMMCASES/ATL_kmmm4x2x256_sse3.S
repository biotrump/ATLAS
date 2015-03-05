/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2015 R. Clint Whaley
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
/*
 * This computational kernel was created using a code fragment demonstrating a
 * Core2Duo-friendly 2-D x86 register block sent to me by Yevgen Voronenko
 * of the CMU/SPIRAL group as a template.  Here is original the code fragment
 * that Yevgen sent me:
 *      movapd    (%rdi), %xmm9
 *      movapd    %xmm9, %xmm6
 *      movapd    48(%rdi), %xmm8
 *      mulpd     %xmm8, %xmm6
 *      addpd     %xmm6, %xmm5
 *      movapd    16(%rdi), %xmm10
 *      movapd    %xmm10, %xmm7
 *      mulpd     %xmm8, %xmm7
 *      addpd     %xmm7, %xmm4
 *      movapd    32(%rdi), %xmm12
 *      mulpd     %xmm12, %xmm8
 *      addpd     %xmm8, %xmm3
 *      movapd    64(%rdi), %xmm11
 *      mulpd     %xmm11, %xmm9
 *      mulpd     %xmm11, %xmm10
 *      mulpd     %xmm11, %xmm12
 *      addpd     %xmm9, %xmm2
 *      addpd     %xmm10, %xmm1
 *      addpd     %xmm12, %xmm0
 */


#if !defined(ATL_GAS_x8664)
   #error "This kernel requires x86-64 assembly!"
#endif

#if !defined(KB) || (KB == 0)
   #error "KB must be a compile-time constant!"
#endif
#if KB > 256
   #error "KB can at most be 256!"
#endif
#if (KB/2)*2 != KB
   #error "KB must be a multiple of 2!"
#endif

#include "atlas_asm.h"
#define movapd movaps
#define nmu     %rdi
#define nnu     %rsi
#define nnu0    %r10
#define pA      %rcx
#define pB      %rax
#define pC      %r9
#define pfA     %rbp
#define pB0     %r12
#define pf0     %rbx
#define pfB     %rdx
#define incAm   %r11
#define pfC     %r13
#define incBn   %r14

#define rA0 	%xmm0
#define rA1 	%xmm1
#define rA2 	%xmm2
#define rA3 	%xmm3
#define rB0 	%xmm4
#define ra0 	%xmm5
#define ra1 	%xmm6
#define ra2 	%xmm7
#define rC00 	%xmm8
#define rC10 	%xmm9
#define rC20 	%xmm10
#define rC30 	%xmm11
#define rC01 	%xmm12
#define rC11 	%xmm13
#define rC21 	%xmm14
#define rC31 	%xmm15

/*
 * Save some inst space by using short version of instructions
 */
#if defined(SREAL) || defined(SCPLX)
   #define movapd movaps
   #define subpd subps
   #define addpd addps
   #define mulpd mulps
   #define addpd addps
   #define kmul 2
   #define SZ 4
#else
   #define movapd movaps
   #define kmul 1
   #define SZ 8
#endif
#ifdef BETAN1
   #define VOP subpd
#elif defined(BETA1)
   #define VOP addpd
#elif defined(VOP)
   #undef VOP
#endif
#ifdef ATL_3DNow
   #define prefC prefetchw
#else
   #define prefC prefetcht0
#endif
#define FSIZE 48
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
      prefetcht0 (pA)
/*
 * Load paramaters
 */
   movq %r8, pB
      prefetcht0 (pB)
   mov nnu, nnu0
   movq FSIZE+16(%rsp), pfB     /* pfB = pBn */
   movq FSIZE+8(%rsp), pfA      /* pfA = pAn */
   movq FSIZE+24(%rsp), pfC     /* pfC = pCn */
      prefetcht0 64(pA)
/*
 * Extend range of small operands by starting at -128
 */
   sub $-128, pA
      prefetcht0 (pA)
   sub $-128, pB
      prefetcht0 64(pA)
   mov $KB*4*SZ, incAm         /* incAm = KB*MU*size */
      prefetcht0 -64(pB)
   mov $KB*2*SZ, incBn         /* incBn = KB*NU*size */
      prefetcht0 (pB)
   movq pB, pB0
   lea (pA, incAm), pf0

   ALIGN16
   .local MNLOOP
   MNLOOP:
/*
 *       Peel first iteration to zero registers
 */
         movapd -128(pA), rC00
         movapd -128(pB), rB0
         movapd rC00, rC01
         mulpd  rB0, rC00
         movapd -112(pA), rC10
         movapd rC10, rC11
         mulpd  rB0, rC10
         movapd -96(pA), rC20
         movapd rC20, rC21
         mulpd  rB0, rC20
         movapd -80(pA), rC30
         movapd rC30, rC31
         mulpd  rB0, rC30
         movapd -112(pB), ra0
         mulpd  ra0, rC01
         mulpd  ra0, rC11
         #if KB > 2*kmul
            movapd -96(pB), rB0
         #endif
         #if KB < 14*kmul
            prefetcht0 (pC)
         #endif
         mulpd  ra0, rC21
         mulpd  ra0, rC31
            prefetcht0 (pf0)
            add $64, pf0
         #if KB > 2*kmul
            movapd -64(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd -48(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd -32(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd -16(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd -80(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 4*kmul
               movapd -64(pB), rB0
            #endif
            #if kmul*(4+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 4*kmul
            movapd 0(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 16(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 32(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 48(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd -48(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 6*kmul
               movapd -32(pB), rB0
            #endif
            #if kmul*(6+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 6*kmul
            movapd 64(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 80(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 96(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 112(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd -16(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 8*kmul
               movapd 0(pB), rB0
            #endif
            #if kmul*(8+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 8*kmul
            movapd 128(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 144(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 160(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 176(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 16(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 10*kmul
               movapd 32(pB), rB0
            #endif
            #if kmul*(10+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 10*kmul
            movapd 192(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 208(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 224(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 240(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 48(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 12*kmul
               movapd 64(pB), rB0
            #endif
            #if kmul*(12+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 12*kmul
            movapd 256(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 272(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 288(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 304(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 80(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 14*kmul
               movapd 96(pB), rB0
            #endif
            #if kmul*(14+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 14*kmul
            movapd 320(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 336(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 352(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 368(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 112(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 16*kmul
               movapd 128(pB), rB0
            #endif
            #if kmul*(16+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 16*kmul
            movapd 384(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 400(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 416(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 432(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 144(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 18*kmul
               movapd 160(pB), rB0
            #endif
            #if kmul*(18+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 18*kmul
            movapd 448(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 464(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 480(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 496(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 176(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 20*kmul
               movapd 192(pB), rB0
            #endif
            #if kmul*(20+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 20*kmul
            movapd 512(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 528(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 544(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 560(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 208(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 22*kmul
               movapd 224(pB), rB0
            #endif
            #if kmul*(22+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 22*kmul
            movapd 576(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 592(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 608(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 624(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 240(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 24*kmul
               movapd 256(pB), rB0
            #endif
            #if kmul*(24+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 24*kmul
            movapd 640(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 656(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 672(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 688(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 272(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 26*kmul
               movapd 288(pB), rB0
            #endif
            #if kmul*(26+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 26*kmul
            movapd 704(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 720(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 736(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 752(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 304(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 28*kmul
               movapd 320(pB), rB0
            #endif
            #if kmul*(28+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 28*kmul
            movapd 768(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 784(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 800(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 816(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 336(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 30*kmul
               movapd 352(pB), rB0
            #endif
            #if kmul*(30+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 30*kmul
            movapd 832(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 848(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 864(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 880(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 368(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 32*kmul
               movapd 384(pB), rB0
            #endif
            #if kmul*(32+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 32*kmul
            movapd 896(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 912(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 928(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 944(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 400(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 34*kmul
               movapd 416(pB), rB0
            #endif
            #if kmul*(34+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 34*kmul
            movapd 960(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 976(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 992(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1008(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 432(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 36*kmul
               movapd 448(pB), rB0
            #endif
            #if kmul*(36+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 36*kmul
            movapd 1024(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1040(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1056(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1072(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 464(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 38*kmul
               movapd 480(pB), rB0
            #endif
            #if kmul*(38+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 38*kmul
            movapd 1088(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1104(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1120(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1136(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 496(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 40*kmul
               movapd 512(pB), rB0
            #endif
            #if kmul*(40+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 40*kmul
            movapd 1152(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1168(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1184(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1200(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 528(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 42*kmul
               movapd 544(pB), rB0
            #endif
            #if kmul*(42+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 42*kmul
            movapd 1216(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1232(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1248(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1264(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 560(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 44*kmul
               movapd 576(pB), rB0
            #endif
            #if kmul*(44+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 44*kmul
            movapd 1280(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1296(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1312(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1328(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 592(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 46*kmul
               movapd 608(pB), rB0
            #endif
            #if kmul*(46+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 46*kmul
            movapd 1344(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1360(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1376(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1392(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 624(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 48*kmul
               movapd 640(pB), rB0
            #endif
            #if kmul*(48+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 48*kmul
            movapd 1408(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1424(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1440(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1456(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 656(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 50*kmul
               movapd 672(pB), rB0
            #endif
            #if kmul*(50+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 50*kmul
            movapd 1472(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1488(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1504(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1520(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 688(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 52*kmul
               movapd 704(pB), rB0
            #endif
            #if kmul*(52+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 52*kmul
            movapd 1536(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1552(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1568(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1584(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 720(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 54*kmul
               movapd 736(pB), rB0
            #endif
            #if kmul*(54+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 54*kmul
            movapd 1600(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1616(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1632(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1648(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 752(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 56*kmul
               movapd 768(pB), rB0
            #endif
            #if kmul*(56+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 56*kmul
            movapd 1664(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1680(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1696(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1712(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 784(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 58*kmul
               movapd 800(pB), rB0
            #endif
            #if kmul*(58+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 58*kmul
            movapd 1728(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1744(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1760(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1776(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 816(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 60*kmul
               movapd 832(pB), rB0
            #endif
            #if kmul*(60+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 60*kmul
            movapd 1792(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1808(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1824(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1840(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 848(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 62*kmul
               movapd 864(pB), rB0
            #endif
            #if kmul*(62+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 62*kmul
            movapd 1856(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1872(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1888(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1904(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 880(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 64*kmul
               movapd 896(pB), rB0
            #endif
            #if kmul*(64+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 64*kmul
            movapd 1920(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 1936(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 1952(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 1968(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 912(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 66*kmul
               movapd 928(pB), rB0
            #endif
            #if kmul*(66+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 66*kmul
            movapd 1984(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2000(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2016(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2032(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 944(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 68*kmul
               movapd 960(pB), rB0
            #endif
            #if kmul*(68+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 68*kmul
            movapd 2048(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2064(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2080(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2096(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 976(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 70*kmul
               movapd 992(pB), rB0
            #endif
            #if kmul*(70+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 70*kmul
            movapd 2112(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2128(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2144(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2160(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1008(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 72*kmul
               movapd 1024(pB), rB0
            #endif
            #if kmul*(72+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 72*kmul
            movapd 2176(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2192(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2208(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2224(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1040(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 74*kmul
               movapd 1056(pB), rB0
            #endif
            #if kmul*(74+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 74*kmul
            movapd 2240(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2256(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2272(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2288(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1072(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 76*kmul
               movapd 1088(pB), rB0
            #endif
            #if kmul*(76+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 76*kmul
            movapd 2304(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2320(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2336(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2352(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1104(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 78*kmul
               movapd 1120(pB), rB0
            #endif
            #if kmul*(78+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 78*kmul
            movapd 2368(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2384(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2400(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2416(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1136(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 80*kmul
               movapd 1152(pB), rB0
            #endif
            #if kmul*(80+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 80*kmul
            movapd 2432(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2448(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2464(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2480(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1168(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 82*kmul
               movapd 1184(pB), rB0
            #endif
            #if kmul*(82+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 82*kmul
            movapd 2496(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2512(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2528(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2544(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1200(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 84*kmul
               movapd 1216(pB), rB0
            #endif
            #if kmul*(84+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 84*kmul
            movapd 2560(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2576(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2592(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2608(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1232(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 86*kmul
               movapd 1248(pB), rB0
            #endif
            #if kmul*(86+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 86*kmul
            movapd 2624(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2640(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2656(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2672(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1264(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 88*kmul
               movapd 1280(pB), rB0
            #endif
            #if kmul*(88+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 88*kmul
            movapd 2688(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2704(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2720(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2736(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1296(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 90*kmul
               movapd 1312(pB), rB0
            #endif
            #if kmul*(90+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 90*kmul
            movapd 2752(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2768(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2784(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2800(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1328(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 92*kmul
               movapd 1344(pB), rB0
            #endif
            #if kmul*(92+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 92*kmul
            movapd 2816(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2832(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2848(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2864(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1360(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 94*kmul
               movapd 1376(pB), rB0
            #endif
            #if kmul*(94+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 94*kmul
            movapd 2880(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2896(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2912(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2928(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1392(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 96*kmul
               movapd 1408(pB), rB0
            #endif
            #if kmul*(96+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 96*kmul
            movapd 2944(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 2960(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 2976(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 2992(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1424(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 98*kmul
               movapd 1440(pB), rB0
            #endif
            #if kmul*(98+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 98*kmul
            movapd 3008(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3024(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3040(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3056(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1456(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 100*kmul
               movapd 1472(pB), rB0
            #endif
            #if kmul*(100+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 100*kmul
            movapd 3072(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3088(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3104(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3120(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1488(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 102*kmul
               movapd 1504(pB), rB0
            #endif
            #if kmul*(102+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 102*kmul
            movapd 3136(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3152(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3168(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3184(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1520(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 104*kmul
               movapd 1536(pB), rB0
            #endif
            #if kmul*(104+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 104*kmul
            movapd 3200(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3216(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3232(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3248(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1552(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 106*kmul
               movapd 1568(pB), rB0
            #endif
            #if kmul*(106+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 106*kmul
            movapd 3264(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3280(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3296(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3312(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1584(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 108*kmul
               movapd 1600(pB), rB0
            #endif
            #if kmul*(108+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 108*kmul
            movapd 3328(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3344(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3360(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3376(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1616(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 110*kmul
               movapd 1632(pB), rB0
            #endif
            #if kmul*(110+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 110*kmul
            movapd 3392(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3408(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3424(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3440(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1648(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 112*kmul
               movapd 1664(pB), rB0
            #endif
            #if kmul*(112+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 112*kmul
            movapd 3456(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3472(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3488(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3504(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1680(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 114*kmul
               movapd 1696(pB), rB0
            #endif
            #if kmul*(114+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 114*kmul
            movapd 3520(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3536(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3552(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3568(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1712(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 116*kmul
               movapd 1728(pB), rB0
            #endif
            #if kmul*(116+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 116*kmul
            movapd 3584(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3600(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3616(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3632(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1744(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 118*kmul
               movapd 1760(pB), rB0
            #endif
            #if kmul*(118+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 118*kmul
            movapd 3648(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3664(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3680(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3696(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1776(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 120*kmul
               movapd 1792(pB), rB0
            #endif
            #if kmul*(120+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 120*kmul
            movapd 3712(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3728(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3744(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3760(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1808(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 122*kmul
               movapd 1824(pB), rB0
            #endif
            #if kmul*(122+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 122*kmul
            movapd 3776(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3792(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3808(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3824(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1840(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 124*kmul
               movapd 1856(pB), rB0
            #endif
            #if kmul*(124+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 124*kmul
            movapd 3840(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3856(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3872(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3888(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1872(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 126*kmul
               movapd 1888(pB), rB0
            #endif
            #if kmul*(126+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 126*kmul
            movapd 3904(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3920(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 3936(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 3952(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1904(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 128*kmul
               movapd 1920(pB), rB0
            #endif
            #if kmul*(128+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 128*kmul
            movapd 3968(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 3984(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4000(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4016(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1936(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 130*kmul
               movapd 1952(pB), rB0
            #endif
            #if kmul*(130+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 130*kmul
            movapd 4032(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4048(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4064(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4080(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 1968(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 132*kmul
               movapd 1984(pB), rB0
            #endif
            #if kmul*(132+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 132*kmul
            movapd 4096(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4112(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4128(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4144(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2000(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 134*kmul
               movapd 2016(pB), rB0
            #endif
            #if kmul*(134+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 134*kmul
            movapd 4160(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4176(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4192(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4208(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2032(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 136*kmul
               movapd 2048(pB), rB0
            #endif
            #if kmul*(136+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 136*kmul
            movapd 4224(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4240(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4256(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4272(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2064(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 138*kmul
               movapd 2080(pB), rB0
            #endif
            #if kmul*(138+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 138*kmul
            movapd 4288(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4304(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4320(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4336(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2096(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 140*kmul
               movapd 2112(pB), rB0
            #endif
            #if kmul*(140+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 140*kmul
            movapd 4352(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4368(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4384(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4400(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2128(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 142*kmul
               movapd 2144(pB), rB0
            #endif
            #if kmul*(142+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 142*kmul
            movapd 4416(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4432(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4448(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4464(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2160(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 144*kmul
               movapd 2176(pB), rB0
            #endif
            #if kmul*(144+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 144*kmul
            movapd 4480(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4496(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4512(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4528(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2192(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 146*kmul
               movapd 2208(pB), rB0
            #endif
            #if kmul*(146+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 146*kmul
            movapd 4544(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4560(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4576(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4592(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2224(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 148*kmul
               movapd 2240(pB), rB0
            #endif
            #if kmul*(148+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 148*kmul
            movapd 4608(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4624(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4640(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4656(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2256(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 150*kmul
               movapd 2272(pB), rB0
            #endif
            #if kmul*(150+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 150*kmul
            movapd 4672(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4688(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4704(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4720(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2288(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 152*kmul
               movapd 2304(pB), rB0
            #endif
            #if kmul*(152+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 152*kmul
            movapd 4736(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4752(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4768(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4784(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2320(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 154*kmul
               movapd 2336(pB), rB0
            #endif
            #if kmul*(154+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 154*kmul
            movapd 4800(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4816(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4832(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4848(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2352(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 156*kmul
               movapd 2368(pB), rB0
            #endif
            #if kmul*(156+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 156*kmul
            movapd 4864(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4880(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4896(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4912(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2384(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 158*kmul
               movapd 2400(pB), rB0
            #endif
            #if kmul*(158+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 158*kmul
            movapd 4928(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 4944(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 4960(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 4976(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2416(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 160*kmul
               movapd 2432(pB), rB0
            #endif
            #if kmul*(160+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 160*kmul
            movapd 4992(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5008(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5024(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5040(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2448(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 162*kmul
               movapd 2464(pB), rB0
            #endif
            #if kmul*(162+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 162*kmul
            movapd 5056(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5072(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5088(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5104(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2480(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 164*kmul
               movapd 2496(pB), rB0
            #endif
            #if kmul*(164+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 164*kmul
            movapd 5120(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5136(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5152(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5168(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2512(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 166*kmul
               movapd 2528(pB), rB0
            #endif
            #if kmul*(166+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 166*kmul
            movapd 5184(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5200(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5216(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5232(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2544(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 168*kmul
               movapd 2560(pB), rB0
            #endif
            #if kmul*(168+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 168*kmul
            movapd 5248(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5264(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5280(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5296(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2576(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 170*kmul
               movapd 2592(pB), rB0
            #endif
            #if kmul*(170+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 170*kmul
            movapd 5312(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5328(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5344(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5360(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2608(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 172*kmul
               movapd 2624(pB), rB0
            #endif
            #if kmul*(172+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 172*kmul
            movapd 5376(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5392(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5408(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5424(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2640(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 174*kmul
               movapd 2656(pB), rB0
            #endif
            #if kmul*(174+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 174*kmul
            movapd 5440(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5456(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5472(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5488(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2672(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 176*kmul
               movapd 2688(pB), rB0
            #endif
            #if kmul*(176+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 176*kmul
            movapd 5504(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5520(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5536(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5552(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2704(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 178*kmul
               movapd 2720(pB), rB0
            #endif
            #if kmul*(178+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 178*kmul
            movapd 5568(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5584(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5600(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5616(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2736(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 180*kmul
               movapd 2752(pB), rB0
            #endif
            #if kmul*(180+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 180*kmul
            movapd 5632(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5648(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5664(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5680(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2768(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 182*kmul
               movapd 2784(pB), rB0
            #endif
            #if kmul*(182+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 182*kmul
            movapd 5696(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5712(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5728(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5744(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2800(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 184*kmul
               movapd 2816(pB), rB0
            #endif
            #if kmul*(184+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 184*kmul
            movapd 5760(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5776(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5792(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5808(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2832(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 186*kmul
               movapd 2848(pB), rB0
            #endif
            #if kmul*(186+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 186*kmul
            movapd 5824(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5840(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5856(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5872(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2864(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 188*kmul
               movapd 2880(pB), rB0
            #endif
            #if kmul*(188+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 188*kmul
            movapd 5888(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5904(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5920(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 5936(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2896(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 190*kmul
               movapd 2912(pB), rB0
            #endif
            #if kmul*(190+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 190*kmul
            movapd 5952(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 5968(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 5984(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6000(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2928(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 192*kmul
               movapd 2944(pB), rB0
            #endif
            #if kmul*(192+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 192*kmul
            movapd 6016(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6032(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6048(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6064(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2960(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 194*kmul
               movapd 2976(pB), rB0
            #endif
            #if kmul*(194+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 194*kmul
            movapd 6080(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6096(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6112(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6128(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 2992(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 196*kmul
               movapd 3008(pB), rB0
            #endif
            #if kmul*(196+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 196*kmul
            movapd 6144(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6160(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6176(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6192(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3024(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 198*kmul
               movapd 3040(pB), rB0
            #endif
            #if kmul*(198+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 198*kmul
            movapd 6208(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6224(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6240(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6256(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3056(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 200*kmul
               movapd 3072(pB), rB0
            #endif
            #if kmul*(200+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 200*kmul
            movapd 6272(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6288(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6304(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6320(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3088(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 202*kmul
               movapd 3104(pB), rB0
            #endif
            #if kmul*(202+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 202*kmul
            movapd 6336(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6352(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6368(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6384(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3120(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 204*kmul
               movapd 3136(pB), rB0
            #endif
            #if kmul*(204+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 204*kmul
            movapd 6400(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6416(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6432(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6448(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3152(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 206*kmul
               movapd 3168(pB), rB0
            #endif
            #if kmul*(206+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 206*kmul
            movapd 6464(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6480(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6496(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6512(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3184(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 208*kmul
               movapd 3200(pB), rB0
            #endif
            #if kmul*(208+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 208*kmul
            movapd 6528(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6544(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6560(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6576(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3216(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 210*kmul
               movapd 3232(pB), rB0
            #endif
            #if kmul*(210+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 210*kmul
            movapd 6592(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6608(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6624(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6640(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3248(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 212*kmul
               movapd 3264(pB), rB0
            #endif
            #if kmul*(212+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 212*kmul
            movapd 6656(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6672(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6688(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6704(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3280(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 214*kmul
               movapd 3296(pB), rB0
            #endif
            #if kmul*(214+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 214*kmul
            movapd 6720(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6736(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6752(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6768(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3312(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 216*kmul
               movapd 3328(pB), rB0
            #endif
            #if kmul*(216+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 216*kmul
            movapd 6784(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6800(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6816(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6832(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3344(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 218*kmul
               movapd 3360(pB), rB0
            #endif
            #if kmul*(218+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 218*kmul
            movapd 6848(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6864(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6880(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6896(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3376(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 220*kmul
               movapd 3392(pB), rB0
            #endif
            #if kmul*(220+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 220*kmul
            movapd 6912(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6928(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 6944(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 6960(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3408(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 222*kmul
               movapd 3424(pB), rB0
            #endif
            #if kmul*(222+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 222*kmul
            movapd 6976(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 6992(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7008(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7024(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3440(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 224*kmul
               movapd 3456(pB), rB0
            #endif
            #if kmul*(224+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 224*kmul
            movapd 7040(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7056(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7072(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7088(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3472(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 226*kmul
               movapd 3488(pB), rB0
            #endif
            #if kmul*(226+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 226*kmul
            movapd 7104(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7120(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7136(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7152(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3504(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 228*kmul
               movapd 3520(pB), rB0
            #endif
            #if kmul*(228+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 228*kmul
            movapd 7168(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7184(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7200(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7216(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3536(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 230*kmul
               movapd 3552(pB), rB0
            #endif
            #if kmul*(230+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 230*kmul
            movapd 7232(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7248(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7264(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7280(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3568(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 232*kmul
               movapd 3584(pB), rB0
            #endif
            #if kmul*(232+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 232*kmul
            movapd 7296(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7312(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7328(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7344(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3600(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 234*kmul
               movapd 3616(pB), rB0
            #endif
            #if kmul*(234+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 234*kmul
            movapd 7360(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7376(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7392(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7408(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3632(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 236*kmul
               movapd 3648(pB), rB0
            #endif
            #if kmul*(236+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 236*kmul
            movapd 7424(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7440(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7456(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7472(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3664(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 238*kmul
               movapd 3680(pB), rB0
            #endif
            #if kmul*(238+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 238*kmul
            movapd 7488(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7504(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7520(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7536(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3696(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 240*kmul
               movapd 3712(pB), rB0
            #endif
            #if kmul*(240+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 240*kmul
            movapd 7552(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7568(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7584(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7600(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3728(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 242*kmul
               movapd 3744(pB), rB0
            #endif
            #if kmul*(242+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 242*kmul
            movapd 7616(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7632(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7648(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7664(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3760(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 244*kmul
               movapd 3776(pB), rB0
            #endif
            #if kmul*(244+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 244*kmul
            movapd 7680(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7696(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7712(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7728(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3792(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 246*kmul
               movapd 3808(pB), rB0
            #endif
            #if kmul*(246+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 246*kmul
            movapd 7744(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7760(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7776(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7792(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3824(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 248*kmul
               movapd 3840(pB), rB0
            #endif
            #if kmul*(248+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 248*kmul
            movapd 7808(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7824(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7840(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7856(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3856(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 250*kmul
               movapd 3872(pB), rB0
            #endif
            #if kmul*(250+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 250*kmul
            movapd 7872(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7888(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7904(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7920(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3888(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 252*kmul
               movapd 3904(pB), rB0
            #endif
            #if kmul*(252+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 252*kmul
            movapd 7936(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 7952(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 7968(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 7984(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3920(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 254*kmul
               movapd 3936(pB), rB0
            #endif
            #if kmul*(254+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
         #if KB > 254*kmul
            movapd 8000(pA), rA0
            movapd rA0, ra0
            mulpd  rB0, rA0
            addpd  rA0, rC00
            movapd 8016(pA), rA1
            movapd rA1, ra1
            mulpd  rB0, rA1
            addpd  rA1, rC10
            movapd 8032(pA), rA2
            movapd rA2, ra2
            mulpd  rB0, rA2
            addpd  rA2, rC20
            movapd 8048(pA), rA3
            mulpd  rA3, rB0
            addpd  rB0, rC30
            movapd 3952(pB), rB0
            mulpd  rB0, ra0
            mulpd  rB0, ra1
            mulpd  rB0, ra2
            mulpd  rB0, rA3
            #if KB > 256*kmul
               movapd 3968(pB), rB0
            #endif
            #if kmul*(256+12) == KB
                prefetcht0 (pC)
            #endif
            addpd  ra0, rC01
            addpd  ra1, rC11
            addpd  ra2, rC21
            addpd  rA3, rC31
         #endif
            prefetcht2 (pfA)
            add $64, pfA
            prefetcht2 (pfB)
            add $64, pfB

/* KDONE: */
/*
 *       Sum up rCx regs, apply original C, and write C out
 */
         #ifdef SREAL
            haddps rC10, rC00  /* rC00 = {c1cd,c1ab,c0cd,c0ab} */
            haddps rC30, rC20  /* rC20 = {c3cd,c3ab,c2cd,c0ab} */
            haddps rC20, rC00  /* rC00 = {c3,c2,c1,c0} */
            #ifdef VOP
               VOP (pC), rC00
            #endif
            movaps rC00, (pC)

            haddps rC11, rC01  /* rC01 = {c1cd,c1ab,c0cd,c0ab} */
            haddps rC31, rC21  /* rC21 = {c3cd,c3ab,c2cd,c0ab} */
            haddps rC21, rC01  /* rC01 = {c3,c2,c1,c0} */
            #ifdef VOP
               VOP 16(pC), rC01
            #endif
            movaps rC01, 16(pC)
            add $32, pC
         #else
            haddpd rC10, rC00  /* rC00 = {c10ab,c0ab} */
            #ifdef VOP
               VOP (pC), rC00
            #endif
            movapd rC00, (pC)
            haddpd rC30, rC20  /* rC20 = {c30ab,c2ab} */
            #ifdef VOP
               VOP 16(pC), rC20
            #endif
            movapd rC20, 16(pC)

            haddpd rC11, rC01
            #ifdef VOP
               VOP 32(pC), rC01
            #endif
            movapd rC01, 32(pC)
            haddpd rC31, rC21
            #ifdef VOP
               VOP 48(pC), rC21
            #endif
            movapd rC21, 48(pC)
            add $64, pC
         #endif
         add incBn, pB   /* KB*vlen*NU*sizeof */
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
   movq    32(%rsp), %r14
   add $FSIZE, %rsp
   ret
