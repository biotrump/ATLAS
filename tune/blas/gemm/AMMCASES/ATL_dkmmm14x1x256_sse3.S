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
#ifndef ATL_GAS_x8664
   #error "This kernel requires x86-64 assembly!"
#endif
#ifndef ATL_SSE3
   #error "This routine requires SSE3!"
#endif
#include "atlas_asm.h"

#define rA0     %xmm0
#define rB0     %xmm1
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
#define rC12    %xmm14
#define rC13    %xmm15

/* #define KK      %rdx */   /* API register */
#define pA      %rcx   /* API reg */
#define pB      %rax   /* comes in as r8 */
#define NMU     %rdi   /* API reg */
#define NNU     %rsi   /* API reg */
#define pC      %r9    /* API reg */
#define pfB     %r10
#define pfA     %r8

#define NNU0    %r11
#define incA    %r12
#define pB0     %r13
#define r112    %r14
/* #define K0      %r15 */
#define FSIZE 6*8
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
#if KB <= 40
   #define prefB(m_) prefetcht0 m_
#else
   #define prefB(m_) prefetcht2 m_
#endif
#if KB > 120
   #define prefA(m_)
#elif KB > 40
   #define prefA(m_) prefetcht2 m_
#else
   #define prefA(m_) prefetcht0 m_
#endif
#define MOVAPD movaps
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

   movq FSIZE+16(%rsp), pfB
   mov %r8, pB
   movq FSIZE+8(%rsp), pfA
   sub $-128, pA
   sub $-128, pB
   mov pB, pB0
   mov $14*KB*8, incA   /* incA = 14*K*sizeof */
   mov NNU, NNU0
   mov $112, r112

   ALIGN16
   MLOOP:
/*      NLOOP: */
         MOVAPD -128(pB), rC13
         MOVAPD -128(pA), rC0
         mulpd rC13, rC0
         MOVAPD -112(pA), rC1
         mulpd rC13, rC1
         MOVAPD -96(pA), rC2
         mulpd rC13, rC2
         MOVAPD -80(pA), rC3
         mulpd rC13, rC3
         MOVAPD -64(pA), rC4
         mulpd rC13, rC4
         MOVAPD -48(pA), rC5
         mulpd rC13, rC5
         MOVAPD -32(pA), rC6
         mulpd rC13, rC6
         MOVAPD -16(pA), rC7
         mulpd rC13, rC7
         MOVAPD (pA), rC8
         mulpd rC13, rC8
         MOVAPD 16(pA), rC9
         mulpd rC13, rC9
         MOVAPD 32(pA), rC10
         mulpd rC13, rC10
         MOVAPD 48(pA), rC11
         mulpd rC13, rC11
         MOVAPD 64(pA), rC12
            prefC (pC)
         mulpd rC13, rC12
            prefC 64(pC)
         mulpd  80(pA), rC13

/*         KLOOP: */
         #if KB > 2
            MOVAPD 96(pA), rA0
            MOVAPD -112(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 288(pA), rA0
            #if KB == 4
               prefB((pfB))
            #elif KB-2 == 4
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 4
               prefB(64(pfB))
            #elif KB-2 == 4
               prefA(64(pfA))
            #endif
            mulpd  304(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 4
            MOVAPD 320(pA), rA0
            MOVAPD -96(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 512(pA), rA0
            #if KB == 6
               prefB((pfB))
            #elif KB-2 == 6
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 6
               prefB(64(pfB))
            #elif KB-2 == 6
               prefA(64(pfA))
            #endif
            mulpd  528(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 6
            MOVAPD 544(pA), rA0
            MOVAPD -80(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 736(pA), rA0
            #if KB == 8
               prefB((pfB))
            #elif KB-2 == 8
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 8
               prefB(64(pfB))
            #elif KB-2 == 8
               prefA(64(pfA))
            #endif
            mulpd  752(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 8
            MOVAPD 768(pA), rA0
            MOVAPD -64(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 960(pA), rA0
            #if KB == 10
               prefB((pfB))
            #elif KB-2 == 10
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 10
               prefB(64(pfB))
            #elif KB-2 == 10
               prefA(64(pfA))
            #endif
            mulpd  976(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 10
            MOVAPD 992(pA), rA0
            MOVAPD -48(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 1008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 1024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 1040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 1056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 1072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 1088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 1104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 1120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 1136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 1152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 1168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 1184(pA), rA0
            #if KB == 12
               prefB((pfB))
            #elif KB-2 == 12
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 12
               prefB(64(pfB))
            #elif KB-2 == 12
               prefA(64(pfA))
            #endif
            mulpd  1200(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 12
            MOVAPD 1216(pA), rA0
            MOVAPD -32(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 1232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 1248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 1264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 1280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 1296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 1312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 1328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 1344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 1360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 1376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 1392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 1408(pA), rA0
            #if KB == 14
               prefB((pfB))
            #elif KB-2 == 14
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 14
               prefB(64(pfB))
            #elif KB-2 == 14
               prefA(64(pfA))
            #endif
            mulpd  1424(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 14
            MOVAPD 1440(pA), rA0
            MOVAPD -16(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 1456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 1472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 1488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 1504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 1520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 1536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 1552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 1568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 1584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 1600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 1616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 1632(pA), rA0
            #if KB == 16
               prefB((pfB))
            #elif KB-2 == 16
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 16
               prefB(64(pfB))
            #elif KB-2 == 16
               prefA(64(pfA))
            #endif
            mulpd  1648(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 16
            MOVAPD 1664(pA), rA0
            MOVAPD 0(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 1680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 1696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 1712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 1728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 1744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 1760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 1776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 1792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 1808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 1824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 1840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 1856(pA), rA0
            #if KB == 18
               prefB((pfB))
            #elif KB-2 == 18
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 18
               prefB(64(pfB))
            #elif KB-2 == 18
               prefA(64(pfA))
            #endif
            mulpd  1872(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 18
            MOVAPD 1888(pA), rA0
            MOVAPD 16(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 1904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 1920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 1936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 1952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 1968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 1984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 2000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 2016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 2032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 2048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 2064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 2080(pA), rA0
            #if KB == 20
               prefB((pfB))
            #elif KB-2 == 20
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 20
               prefB(64(pfB))
            #elif KB-2 == 20
               prefA(64(pfA))
            #endif
            mulpd  2096(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 20
            MOVAPD 2112(pA), rA0
            MOVAPD 32(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 2128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 2144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 2160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 2176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 2192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 2208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 2224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 2240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 2256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 2272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 2288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 2304(pA), rA0
            #if KB == 22
               prefB((pfB))
            #elif KB-2 == 22
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 22
               prefB(64(pfB))
            #elif KB-2 == 22
               prefA(64(pfA))
            #endif
            mulpd  2320(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 22
            MOVAPD 2336(pA), rA0
            MOVAPD 48(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 2352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 2368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 2384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 2400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 2416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 2432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 2448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 2464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 2480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 2496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 2512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 2528(pA), rA0
            #if KB == 24
               prefB((pfB))
            #elif KB-2 == 24
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 24
               prefB(64(pfB))
            #elif KB-2 == 24
               prefA(64(pfA))
            #endif
            mulpd  2544(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 24
            MOVAPD 2560(pA), rA0
            MOVAPD 64(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 2576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 2592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 2608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 2624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 2640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 2656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 2672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 2688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 2704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 2720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 2736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 2752(pA), rA0
            #if KB == 26
               prefB((pfB))
            #elif KB-2 == 26
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 26
               prefB(64(pfB))
            #elif KB-2 == 26
               prefA(64(pfA))
            #endif
            mulpd  2768(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 26
            MOVAPD 2784(pA), rA0
            MOVAPD 80(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 2800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 2816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 2832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 2848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 2864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 2880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 2896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 2912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 2928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 2944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 2960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 2976(pA), rA0
            #if KB == 28
               prefB((pfB))
            #elif KB-2 == 28
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 28
               prefB(64(pfB))
            #elif KB-2 == 28
               prefA(64(pfA))
            #endif
            mulpd  2992(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 28
            MOVAPD 3008(pA), rA0
            MOVAPD 96(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 3024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 3040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 3056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 3072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 3088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 3104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 3120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 3136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 3152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 3168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 3184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 3200(pA), rA0
            #if KB == 30
               prefB((pfB))
            #elif KB-2 == 30
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 30
               prefB(64(pfB))
            #elif KB-2 == 30
               prefA(64(pfA))
            #endif
            mulpd  3216(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 30
            MOVAPD 3232(pA), rA0
            MOVAPD 112(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 3248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 3264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 3280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 3296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 3312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 3328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 3344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 3360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 3376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 3392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 3408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 3424(pA), rA0
            #if KB == 32
               prefB((pfB))
            #elif KB-2 == 32
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 32
               prefB(64(pfB))
            #elif KB-2 == 32
               prefA(64(pfA))
            #endif
            mulpd  3440(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 32
            MOVAPD 3456(pA), rA0
            MOVAPD 128(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 3472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 3488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 3504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 3520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 3536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 3552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 3568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 3584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 3600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 3616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 3632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 3648(pA), rA0
            #if KB == 34
               prefB((pfB))
            #elif KB-2 == 34
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 34
               prefB(64(pfB))
            #elif KB-2 == 34
               prefA(64(pfA))
            #endif
            mulpd  3664(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 34
            MOVAPD 3680(pA), rA0
            MOVAPD 144(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 3696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 3712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 3728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 3744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 3760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 3776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 3792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 3808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 3824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 3840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 3856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 3872(pA), rA0
            #if KB == 36
               prefB((pfB))
            #elif KB-2 == 36
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 36
               prefB(64(pfB))
            #elif KB-2 == 36
               prefA(64(pfA))
            #endif
            mulpd  3888(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 36
            MOVAPD 3904(pA), rA0
            MOVAPD 160(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 3920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 3936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 3952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 3968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 3984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 4000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 4016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 4032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 4048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 4064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 4080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 4096(pA), rA0
            #if KB == 38
               prefB((pfB))
            #elif KB-2 == 38
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 38
               prefB(64(pfB))
            #elif KB-2 == 38
               prefA(64(pfA))
            #endif
            mulpd  4112(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 38
            MOVAPD 4128(pA), rA0
            MOVAPD 176(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 4144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 4160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 4176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 4192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 4208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 4224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 4240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 4256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 4272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 4288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 4304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 4320(pA), rA0
            #if KB == 40
               prefB((pfB))
            #elif KB-2 == 40
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 40
               prefB(64(pfB))
            #elif KB-2 == 40
               prefA(64(pfA))
            #endif
            mulpd  4336(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 40
            MOVAPD 4352(pA), rA0
            MOVAPD 192(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 4368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 4384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 4400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 4416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 4432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 4448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 4464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 4480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 4496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 4512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 4528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 4544(pA), rA0
            #if KB == 42
               prefB((pfB))
            #elif KB-2 == 42
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 42
               prefB(64(pfB))
            #elif KB-2 == 42
               prefA(64(pfA))
            #endif
            mulpd  4560(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 42
            MOVAPD 4576(pA), rA0
            MOVAPD 208(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 4592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 4608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 4624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 4640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 4656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 4672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 4688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 4704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 4720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 4736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 4752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 4768(pA), rA0
            #if KB == 44
               prefB((pfB))
            #elif KB-2 == 44
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 44
               prefB(64(pfB))
            #elif KB-2 == 44
               prefA(64(pfA))
            #endif
            mulpd  4784(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 44
            MOVAPD 4800(pA), rA0
            MOVAPD 224(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 4816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 4832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 4848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 4864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 4880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 4896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 4912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 4928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 4944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 4960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 4976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 4992(pA), rA0
            #if KB == 46
               prefB((pfB))
            #elif KB-2 == 46
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 46
               prefB(64(pfB))
            #elif KB-2 == 46
               prefA(64(pfA))
            #endif
            mulpd  5008(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 46
            MOVAPD 5024(pA), rA0
            MOVAPD 240(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 5040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 5056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 5072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 5088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 5104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 5120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 5136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 5152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 5168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 5184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 5200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 5216(pA), rA0
            #if KB == 48
               prefB((pfB))
            #elif KB-2 == 48
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 48
               prefB(64(pfB))
            #elif KB-2 == 48
               prefA(64(pfA))
            #endif
            mulpd  5232(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 48
            MOVAPD 5248(pA), rA0
            MOVAPD 256(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 5264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 5280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 5296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 5312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 5328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 5344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 5360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 5376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 5392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 5408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 5424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 5440(pA), rA0
            #if KB == 50
               prefB((pfB))
            #elif KB-2 == 50
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 50
               prefB(64(pfB))
            #elif KB-2 == 50
               prefA(64(pfA))
            #endif
            mulpd  5456(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 50
            MOVAPD 5472(pA), rA0
            MOVAPD 272(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 5488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 5504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 5520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 5536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 5552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 5568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 5584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 5600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 5616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 5632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 5648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 5664(pA), rA0
            #if KB == 52
               prefB((pfB))
            #elif KB-2 == 52
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 52
               prefB(64(pfB))
            #elif KB-2 == 52
               prefA(64(pfA))
            #endif
            mulpd  5680(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 52
            MOVAPD 5696(pA), rA0
            MOVAPD 288(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 5712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 5728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 5744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 5760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 5776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 5792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 5808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 5824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 5840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 5856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 5872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 5888(pA), rA0
            #if KB == 54
               prefB((pfB))
            #elif KB-2 == 54
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 54
               prefB(64(pfB))
            #elif KB-2 == 54
               prefA(64(pfA))
            #endif
            mulpd  5904(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 54
            MOVAPD 5920(pA), rA0
            MOVAPD 304(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 5936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 5952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 5968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 5984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 6000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 6016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 6032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 6048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 6064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 6080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 6096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 6112(pA), rA0
            #if KB == 56
               prefB((pfB))
            #elif KB-2 == 56
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 56
               prefB(64(pfB))
            #elif KB-2 == 56
               prefA(64(pfA))
            #endif
            mulpd  6128(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 56
            MOVAPD 6144(pA), rA0
            MOVAPD 320(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 6160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 6176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 6192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 6208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 6224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 6240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 6256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 6272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 6288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 6304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 6320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 6336(pA), rA0
            #if KB == 58
               prefB((pfB))
            #elif KB-2 == 58
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 58
               prefB(64(pfB))
            #elif KB-2 == 58
               prefA(64(pfA))
            #endif
            mulpd  6352(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 58
            MOVAPD 6368(pA), rA0
            MOVAPD 336(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 6384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 6400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 6416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 6432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 6448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 6464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 6480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 6496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 6512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 6528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 6544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 6560(pA), rA0
            #if KB == 60
               prefB((pfB))
            #elif KB-2 == 60
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 60
               prefB(64(pfB))
            #elif KB-2 == 60
               prefA(64(pfA))
            #endif
            mulpd  6576(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 60
            MOVAPD 6592(pA), rA0
            MOVAPD 352(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 6608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 6624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 6640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 6656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 6672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 6688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 6704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 6720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 6736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 6752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 6768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 6784(pA), rA0
            #if KB == 62
               prefB((pfB))
            #elif KB-2 == 62
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 62
               prefB(64(pfB))
            #elif KB-2 == 62
               prefA(64(pfA))
            #endif
            mulpd  6800(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 62
            MOVAPD 6816(pA), rA0
            MOVAPD 368(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 6832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 6848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 6864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 6880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 6896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 6912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 6928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 6944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 6960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 6976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 6992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 7008(pA), rA0
            #if KB == 64
               prefB((pfB))
            #elif KB-2 == 64
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 64
               prefB(64(pfB))
            #elif KB-2 == 64
               prefA(64(pfA))
            #endif
            mulpd  7024(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 64
            MOVAPD 7040(pA), rA0
            MOVAPD 384(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 7056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 7072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 7088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 7104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 7120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 7136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 7152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 7168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 7184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 7200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 7216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 7232(pA), rA0
            #if KB == 66
               prefB((pfB))
            #elif KB-2 == 66
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 66
               prefB(64(pfB))
            #elif KB-2 == 66
               prefA(64(pfA))
            #endif
            mulpd  7248(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 66
            MOVAPD 7264(pA), rA0
            MOVAPD 400(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 7280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 7296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 7312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 7328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 7344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 7360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 7376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 7392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 7408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 7424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 7440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 7456(pA), rA0
            #if KB == 68
               prefB((pfB))
            #elif KB-2 == 68
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 68
               prefB(64(pfB))
            #elif KB-2 == 68
               prefA(64(pfA))
            #endif
            mulpd  7472(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 68
            MOVAPD 7488(pA), rA0
            MOVAPD 416(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 7504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 7520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 7536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 7552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 7568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 7584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 7600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 7616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 7632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 7648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 7664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 7680(pA), rA0
            #if KB == 70
               prefB((pfB))
            #elif KB-2 == 70
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 70
               prefB(64(pfB))
            #elif KB-2 == 70
               prefA(64(pfA))
            #endif
            mulpd  7696(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 70
            MOVAPD 7712(pA), rA0
            MOVAPD 432(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 7728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 7744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 7760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 7776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 7792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 7808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 7824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 7840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 7856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 7872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 7888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 7904(pA), rA0
            #if KB == 72
               prefB((pfB))
            #elif KB-2 == 72
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 72
               prefB(64(pfB))
            #elif KB-2 == 72
               prefA(64(pfA))
            #endif
            mulpd  7920(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 72
            MOVAPD 7936(pA), rA0
            MOVAPD 448(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 7952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 7968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 7984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 8000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 8016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 8032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 8048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 8064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 8080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 8096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 8112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 8128(pA), rA0
            #if KB == 74
               prefB((pfB))
            #elif KB-2 == 74
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 74
               prefB(64(pfB))
            #elif KB-2 == 74
               prefA(64(pfA))
            #endif
            mulpd  8144(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 74
            MOVAPD 8160(pA), rA0
            MOVAPD 464(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 8176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 8192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 8208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 8224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 8240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 8256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 8272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 8288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 8304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 8320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 8336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 8352(pA), rA0
            #if KB == 76
               prefB((pfB))
            #elif KB-2 == 76
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 76
               prefB(64(pfB))
            #elif KB-2 == 76
               prefA(64(pfA))
            #endif
            mulpd  8368(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 76
            MOVAPD 8384(pA), rA0
            MOVAPD 480(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 8400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 8416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 8432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 8448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 8464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 8480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 8496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 8512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 8528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 8544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 8560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 8576(pA), rA0
            #if KB == 78
               prefB((pfB))
            #elif KB-2 == 78
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 78
               prefB(64(pfB))
            #elif KB-2 == 78
               prefA(64(pfA))
            #endif
            mulpd  8592(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 78
            MOVAPD 8608(pA), rA0
            MOVAPD 496(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 8624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 8640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 8656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 8672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 8688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 8704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 8720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 8736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 8752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 8768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 8784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 8800(pA), rA0
            #if KB == 80
               prefB((pfB))
            #elif KB-2 == 80
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 80
               prefB(64(pfB))
            #elif KB-2 == 80
               prefA(64(pfA))
            #endif
            mulpd  8816(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 80
            MOVAPD 8832(pA), rA0
            MOVAPD 512(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 8848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 8864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 8880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 8896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 8912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 8928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 8944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 8960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 8976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 8992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 9008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 9024(pA), rA0
            #if KB == 82
               prefB((pfB))
            #elif KB-2 == 82
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 82
               prefB(64(pfB))
            #elif KB-2 == 82
               prefA(64(pfA))
            #endif
            mulpd  9040(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 82
            MOVAPD 9056(pA), rA0
            MOVAPD 528(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 9072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 9088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 9104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 9120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 9136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 9152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 9168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 9184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 9200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 9216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 9232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 9248(pA), rA0
            #if KB == 84
               prefB((pfB))
            #elif KB-2 == 84
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 84
               prefB(64(pfB))
            #elif KB-2 == 84
               prefA(64(pfA))
            #endif
            mulpd  9264(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 84
            MOVAPD 9280(pA), rA0
            MOVAPD 544(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 9296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 9312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 9328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 9344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 9360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 9376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 9392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 9408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 9424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 9440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 9456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 9472(pA), rA0
            #if KB == 86
               prefB((pfB))
            #elif KB-2 == 86
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 86
               prefB(64(pfB))
            #elif KB-2 == 86
               prefA(64(pfA))
            #endif
            mulpd  9488(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 86
            MOVAPD 9504(pA), rA0
            MOVAPD 560(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 9520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 9536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 9552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 9568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 9584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 9600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 9616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 9632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 9648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 9664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 9680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 9696(pA), rA0
            #if KB == 88
               prefB((pfB))
            #elif KB-2 == 88
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 88
               prefB(64(pfB))
            #elif KB-2 == 88
               prefA(64(pfA))
            #endif
            mulpd  9712(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 88
            MOVAPD 9728(pA), rA0
            MOVAPD 576(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 9744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 9760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 9776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 9792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 9808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 9824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 9840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 9856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 9872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 9888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 9904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 9920(pA), rA0
            #if KB == 90
               prefB((pfB))
            #elif KB-2 == 90
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 90
               prefB(64(pfB))
            #elif KB-2 == 90
               prefA(64(pfA))
            #endif
            mulpd  9936(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 90
            MOVAPD 9952(pA), rA0
            MOVAPD 592(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 9968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 9984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 10000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 10016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 10032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 10048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 10064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 10080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 10096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 10112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 10128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 10144(pA), rA0
            #if KB == 92
               prefB((pfB))
            #elif KB-2 == 92
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 92
               prefB(64(pfB))
            #elif KB-2 == 92
               prefA(64(pfA))
            #endif
            mulpd  10160(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 92
            MOVAPD 10176(pA), rA0
            MOVAPD 608(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 10192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 10208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 10224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 10240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 10256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 10272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 10288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 10304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 10320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 10336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 10352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 10368(pA), rA0
            #if KB == 94
               prefB((pfB))
            #elif KB-2 == 94
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 94
               prefB(64(pfB))
            #elif KB-2 == 94
               prefA(64(pfA))
            #endif
            mulpd  10384(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 94
            MOVAPD 10400(pA), rA0
            MOVAPD 624(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 10416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 10432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 10448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 10464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 10480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 10496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 10512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 10528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 10544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 10560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 10576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 10592(pA), rA0
            #if KB == 96
               prefB((pfB))
            #elif KB-2 == 96
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 96
               prefB(64(pfB))
            #elif KB-2 == 96
               prefA(64(pfA))
            #endif
            mulpd  10608(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 96
            MOVAPD 10624(pA), rA0
            MOVAPD 640(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 10640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 10656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 10672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 10688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 10704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 10720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 10736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 10752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 10768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 10784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 10800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 10816(pA), rA0
            #if KB == 98
               prefB((pfB))
            #elif KB-2 == 98
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 98
               prefB(64(pfB))
            #elif KB-2 == 98
               prefA(64(pfA))
            #endif
            mulpd  10832(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 98
            MOVAPD 10848(pA), rA0
            MOVAPD 656(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 10864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 10880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 10896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 10912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 10928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 10944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 10960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 10976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 10992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 11008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 11024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 11040(pA), rA0
            #if KB == 100
               prefB((pfB))
            #elif KB-2 == 100
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 100
               prefB(64(pfB))
            #elif KB-2 == 100
               prefA(64(pfA))
            #endif
            mulpd  11056(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 100
            MOVAPD 11072(pA), rA0
            MOVAPD 672(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 11088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 11104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 11120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 11136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 11152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 11168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 11184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 11200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 11216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 11232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 11248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 11264(pA), rA0
            #if KB == 102
               prefB((pfB))
            #elif KB-2 == 102
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 102
               prefB(64(pfB))
            #elif KB-2 == 102
               prefA(64(pfA))
            #endif
            mulpd  11280(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 102
            MOVAPD 11296(pA), rA0
            MOVAPD 688(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 11312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 11328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 11344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 11360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 11376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 11392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 11408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 11424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 11440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 11456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 11472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 11488(pA), rA0
            #if KB == 104
               prefB((pfB))
            #elif KB-2 == 104
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 104
               prefB(64(pfB))
            #elif KB-2 == 104
               prefA(64(pfA))
            #endif
            mulpd  11504(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 104
            MOVAPD 11520(pA), rA0
            MOVAPD 704(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 11536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 11552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 11568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 11584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 11600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 11616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 11632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 11648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 11664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 11680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 11696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 11712(pA), rA0
            #if KB == 106
               prefB((pfB))
            #elif KB-2 == 106
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 106
               prefB(64(pfB))
            #elif KB-2 == 106
               prefA(64(pfA))
            #endif
            mulpd  11728(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 106
            MOVAPD 11744(pA), rA0
            MOVAPD 720(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 11760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 11776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 11792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 11808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 11824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 11840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 11856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 11872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 11888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 11904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 11920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 11936(pA), rA0
            #if KB == 108
               prefB((pfB))
            #elif KB-2 == 108
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 108
               prefB(64(pfB))
            #elif KB-2 == 108
               prefA(64(pfA))
            #endif
            mulpd  11952(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 108
            MOVAPD 11968(pA), rA0
            MOVAPD 736(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 11984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 12000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 12016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 12032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 12048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 12064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 12080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 12096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 12112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 12128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 12144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 12160(pA), rA0
            #if KB == 110
               prefB((pfB))
            #elif KB-2 == 110
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 110
               prefB(64(pfB))
            #elif KB-2 == 110
               prefA(64(pfA))
            #endif
            mulpd  12176(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 110
            MOVAPD 12192(pA), rA0
            MOVAPD 752(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 12208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 12224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 12240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 12256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 12272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 12288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 12304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 12320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 12336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 12352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 12368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 12384(pA), rA0
            #if KB == 112
               prefB((pfB))
            #elif KB-2 == 112
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 112
               prefB(64(pfB))
            #elif KB-2 == 112
               prefA(64(pfA))
            #endif
            mulpd  12400(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 112
            MOVAPD 12416(pA), rA0
            MOVAPD 768(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 12432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 12448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 12464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 12480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 12496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 12512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 12528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 12544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 12560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 12576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 12592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 12608(pA), rA0
            #if KB == 114
               prefB((pfB))
            #elif KB-2 == 114
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 114
               prefB(64(pfB))
            #elif KB-2 == 114
               prefA(64(pfA))
            #endif
            mulpd  12624(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 114
            MOVAPD 12640(pA), rA0
            MOVAPD 784(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 12656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 12672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 12688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 12704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 12720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 12736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 12752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 12768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 12784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 12800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 12816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 12832(pA), rA0
            #if KB == 116
               prefB((pfB))
            #elif KB-2 == 116
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 116
               prefB(64(pfB))
            #elif KB-2 == 116
               prefA(64(pfA))
            #endif
            mulpd  12848(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 116
            MOVAPD 12864(pA), rA0
            MOVAPD 800(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 12880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 12896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 12912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 12928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 12944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 12960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 12976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 12992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 13008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 13024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 13040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 13056(pA), rA0
            #if KB == 118
               prefB((pfB))
            #elif KB-2 == 118
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 118
               prefB(64(pfB))
            #elif KB-2 == 118
               prefA(64(pfA))
            #endif
            mulpd  13072(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 118
            MOVAPD 13088(pA), rA0
            MOVAPD 816(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 13104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 13120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 13136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 13152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 13168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 13184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 13200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 13216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 13232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 13248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 13264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 13280(pA), rA0
            #if KB == 120
               prefB((pfB))
            #elif KB-2 == 120
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 120
               prefB(64(pfB))
            #elif KB-2 == 120
               prefA(64(pfA))
            #endif
            mulpd  13296(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 120
            MOVAPD 13312(pA), rA0
            MOVAPD 832(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 13328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 13344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 13360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 13376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 13392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 13408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 13424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 13440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 13456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 13472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 13488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 13504(pA), rA0
            #if KB == 122
               prefB((pfB))
            #elif KB-2 == 122
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 122
               prefB(64(pfB))
            #elif KB-2 == 122
               prefA(64(pfA))
            #endif
            mulpd  13520(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 122
            MOVAPD 13536(pA), rA0
            MOVAPD 848(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 13552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 13568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 13584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 13600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 13616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 13632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 13648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 13664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 13680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 13696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 13712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 13728(pA), rA0
            #if KB == 124
               prefB((pfB))
            #elif KB-2 == 124
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 124
               prefB(64(pfB))
            #elif KB-2 == 124
               prefA(64(pfA))
            #endif
            mulpd  13744(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 124
            MOVAPD 13760(pA), rA0
            MOVAPD 864(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 13776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 13792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 13808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 13824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 13840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 13856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 13872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 13888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 13904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 13920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 13936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 13952(pA), rA0
            #if KB == 126
               prefB((pfB))
            #elif KB-2 == 126
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 126
               prefB(64(pfB))
            #elif KB-2 == 126
               prefA(64(pfA))
            #endif
            mulpd  13968(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 126
            MOVAPD 13984(pA), rA0
            MOVAPD 880(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 14000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 14016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 14032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 14048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 14064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 14080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 14096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 14112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 14128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 14144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 14160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 14176(pA), rA0
            #if KB == 128
               prefB((pfB))
            #elif KB-2 == 128
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 128
               prefB(64(pfB))
            #elif KB-2 == 128
               prefA(64(pfA))
            #endif
            mulpd  14192(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 128
            MOVAPD 14208(pA), rA0
            MOVAPD 896(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 14224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 14240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 14256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 14272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 14288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 14304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 14320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 14336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 14352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 14368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 14384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 14400(pA), rA0
            #if KB == 130
               prefB((pfB))
            #elif KB-2 == 130
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 130
               prefB(64(pfB))
            #elif KB-2 == 130
               prefA(64(pfA))
            #endif
            mulpd  14416(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 130
            MOVAPD 14432(pA), rA0
            MOVAPD 912(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 14448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 14464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 14480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 14496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 14512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 14528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 14544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 14560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 14576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 14592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 14608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 14624(pA), rA0
            #if KB == 132
               prefB((pfB))
            #elif KB-2 == 132
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 132
               prefB(64(pfB))
            #elif KB-2 == 132
               prefA(64(pfA))
            #endif
            mulpd  14640(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 132
            MOVAPD 14656(pA), rA0
            MOVAPD 928(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 14672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 14688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 14704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 14720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 14736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 14752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 14768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 14784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 14800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 14816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 14832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 14848(pA), rA0
            #if KB == 134
               prefB((pfB))
            #elif KB-2 == 134
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 134
               prefB(64(pfB))
            #elif KB-2 == 134
               prefA(64(pfA))
            #endif
            mulpd  14864(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 134
            MOVAPD 14880(pA), rA0
            MOVAPD 944(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 14896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 14912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 14928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 14944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 14960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 14976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 14992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 15008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 15024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 15040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 15056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 15072(pA), rA0
            #if KB == 136
               prefB((pfB))
            #elif KB-2 == 136
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 136
               prefB(64(pfB))
            #elif KB-2 == 136
               prefA(64(pfA))
            #endif
            mulpd  15088(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 136
            MOVAPD 15104(pA), rA0
            MOVAPD 960(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 15120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 15136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 15152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 15168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 15184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 15200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 15216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 15232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 15248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 15264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 15280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 15296(pA), rA0
            #if KB == 138
               prefB((pfB))
            #elif KB-2 == 138
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 138
               prefB(64(pfB))
            #elif KB-2 == 138
               prefA(64(pfA))
            #endif
            mulpd  15312(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 138
            MOVAPD 15328(pA), rA0
            MOVAPD 976(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 15344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 15360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 15376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 15392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 15408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 15424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 15440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 15456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 15472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 15488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 15504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 15520(pA), rA0
            #if KB == 140
               prefB((pfB))
            #elif KB-2 == 140
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 140
               prefB(64(pfB))
            #elif KB-2 == 140
               prefA(64(pfA))
            #endif
            mulpd  15536(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 140
            MOVAPD 15552(pA), rA0
            MOVAPD 992(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 15568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 15584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 15600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 15616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 15632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 15648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 15664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 15680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 15696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 15712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 15728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 15744(pA), rA0
            #if KB == 142
               prefB((pfB))
            #elif KB-2 == 142
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 142
               prefB(64(pfB))
            #elif KB-2 == 142
               prefA(64(pfA))
            #endif
            mulpd  15760(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 142
            MOVAPD 15776(pA), rA0
            MOVAPD 1008(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 15792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 15808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 15824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 15840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 15856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 15872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 15888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 15904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 15920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 15936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 15952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 15968(pA), rA0
            #if KB == 144
               prefB((pfB))
            #elif KB-2 == 144
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 144
               prefB(64(pfB))
            #elif KB-2 == 144
               prefA(64(pfA))
            #endif
            mulpd  15984(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 144
            MOVAPD 16000(pA), rA0
            MOVAPD 1024(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 16016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 16032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 16048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 16064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 16080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 16096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 16112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 16128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 16144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 16160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 16176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 16192(pA), rA0
            #if KB == 146
               prefB((pfB))
            #elif KB-2 == 146
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 146
               prefB(64(pfB))
            #elif KB-2 == 146
               prefA(64(pfA))
            #endif
            mulpd  16208(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 146
            MOVAPD 16224(pA), rA0
            MOVAPD 1040(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 16240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 16256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 16272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 16288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 16304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 16320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 16336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 16352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 16368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 16384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 16400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 16416(pA), rA0
            #if KB == 148
               prefB((pfB))
            #elif KB-2 == 148
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 148
               prefB(64(pfB))
            #elif KB-2 == 148
               prefA(64(pfA))
            #endif
            mulpd  16432(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 148
            MOVAPD 16448(pA), rA0
            MOVAPD 1056(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 16464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 16480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 16496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 16512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 16528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 16544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 16560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 16576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 16592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 16608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 16624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 16640(pA), rA0
            #if KB == 150
               prefB((pfB))
            #elif KB-2 == 150
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 150
               prefB(64(pfB))
            #elif KB-2 == 150
               prefA(64(pfA))
            #endif
            mulpd  16656(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 150
            MOVAPD 16672(pA), rA0
            MOVAPD 1072(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 16688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 16704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 16720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 16736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 16752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 16768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 16784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 16800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 16816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 16832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 16848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 16864(pA), rA0
            #if KB == 152
               prefB((pfB))
            #elif KB-2 == 152
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 152
               prefB(64(pfB))
            #elif KB-2 == 152
               prefA(64(pfA))
            #endif
            mulpd  16880(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 152
            MOVAPD 16896(pA), rA0
            MOVAPD 1088(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 16912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 16928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 16944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 16960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 16976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 16992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 17008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 17024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 17040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 17056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 17072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 17088(pA), rA0
            #if KB == 154
               prefB((pfB))
            #elif KB-2 == 154
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 154
               prefB(64(pfB))
            #elif KB-2 == 154
               prefA(64(pfA))
            #endif
            mulpd  17104(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 154
            MOVAPD 17120(pA), rA0
            MOVAPD 1104(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 17136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 17152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 17168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 17184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 17200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 17216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 17232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 17248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 17264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 17280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 17296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 17312(pA), rA0
            #if KB == 156
               prefB((pfB))
            #elif KB-2 == 156
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 156
               prefB(64(pfB))
            #elif KB-2 == 156
               prefA(64(pfA))
            #endif
            mulpd  17328(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 156
            MOVAPD 17344(pA), rA0
            MOVAPD 1120(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 17360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 17376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 17392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 17408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 17424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 17440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 17456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 17472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 17488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 17504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 17520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 17536(pA), rA0
            #if KB == 158
               prefB((pfB))
            #elif KB-2 == 158
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 158
               prefB(64(pfB))
            #elif KB-2 == 158
               prefA(64(pfA))
            #endif
            mulpd  17552(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 158
            MOVAPD 17568(pA), rA0
            MOVAPD 1136(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 17584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 17600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 17616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 17632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 17648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 17664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 17680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 17696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 17712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 17728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 17744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 17760(pA), rA0
            #if KB == 160
               prefB((pfB))
            #elif KB-2 == 160
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 160
               prefB(64(pfB))
            #elif KB-2 == 160
               prefA(64(pfA))
            #endif
            mulpd  17776(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 160
            MOVAPD 17792(pA), rA0
            MOVAPD 1152(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 17808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 17824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 17840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 17856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 17872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 17888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 17904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 17920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 17936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 17952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 17968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 17984(pA), rA0
            #if KB == 162
               prefB((pfB))
            #elif KB-2 == 162
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 162
               prefB(64(pfB))
            #elif KB-2 == 162
               prefA(64(pfA))
            #endif
            mulpd  18000(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 162
            MOVAPD 18016(pA), rA0
            MOVAPD 1168(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 18032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 18048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 18064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 18080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 18096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 18112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 18128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 18144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 18160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 18176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 18192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 18208(pA), rA0
            #if KB == 164
               prefB((pfB))
            #elif KB-2 == 164
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 164
               prefB(64(pfB))
            #elif KB-2 == 164
               prefA(64(pfA))
            #endif
            mulpd  18224(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 164
            MOVAPD 18240(pA), rA0
            MOVAPD 1184(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 18256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 18272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 18288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 18304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 18320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 18336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 18352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 18368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 18384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 18400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 18416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 18432(pA), rA0
            #if KB == 166
               prefB((pfB))
            #elif KB-2 == 166
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 166
               prefB(64(pfB))
            #elif KB-2 == 166
               prefA(64(pfA))
            #endif
            mulpd  18448(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 166
            MOVAPD 18464(pA), rA0
            MOVAPD 1200(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 18480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 18496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 18512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 18528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 18544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 18560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 18576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 18592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 18608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 18624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 18640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 18656(pA), rA0
            #if KB == 168
               prefB((pfB))
            #elif KB-2 == 168
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 168
               prefB(64(pfB))
            #elif KB-2 == 168
               prefA(64(pfA))
            #endif
            mulpd  18672(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 168
            MOVAPD 18688(pA), rA0
            MOVAPD 1216(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 18704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 18720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 18736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 18752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 18768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 18784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 18800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 18816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 18832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 18848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 18864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 18880(pA), rA0
            #if KB == 170
               prefB((pfB))
            #elif KB-2 == 170
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 170
               prefB(64(pfB))
            #elif KB-2 == 170
               prefA(64(pfA))
            #endif
            mulpd  18896(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 170
            MOVAPD 18912(pA), rA0
            MOVAPD 1232(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 18928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 18944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 18960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 18976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 18992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 19008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 19024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 19040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 19056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 19072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 19088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 19104(pA), rA0
            #if KB == 172
               prefB((pfB))
            #elif KB-2 == 172
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 172
               prefB(64(pfB))
            #elif KB-2 == 172
               prefA(64(pfA))
            #endif
            mulpd  19120(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 172
            MOVAPD 19136(pA), rA0
            MOVAPD 1248(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 19152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 19168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 19184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 19200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 19216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 19232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 19248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 19264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 19280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 19296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 19312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 19328(pA), rA0
            #if KB == 174
               prefB((pfB))
            #elif KB-2 == 174
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 174
               prefB(64(pfB))
            #elif KB-2 == 174
               prefA(64(pfA))
            #endif
            mulpd  19344(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 174
            MOVAPD 19360(pA), rA0
            MOVAPD 1264(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 19376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 19392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 19408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 19424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 19440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 19456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 19472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 19488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 19504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 19520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 19536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 19552(pA), rA0
            #if KB == 176
               prefB((pfB))
            #elif KB-2 == 176
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 176
               prefB(64(pfB))
            #elif KB-2 == 176
               prefA(64(pfA))
            #endif
            mulpd  19568(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 176
            MOVAPD 19584(pA), rA0
            MOVAPD 1280(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 19600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 19616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 19632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 19648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 19664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 19680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 19696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 19712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 19728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 19744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 19760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 19776(pA), rA0
            #if KB == 178
               prefB((pfB))
            #elif KB-2 == 178
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 178
               prefB(64(pfB))
            #elif KB-2 == 178
               prefA(64(pfA))
            #endif
            mulpd  19792(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 178
            MOVAPD 19808(pA), rA0
            MOVAPD 1296(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 19824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 19840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 19856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 19872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 19888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 19904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 19920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 19936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 19952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 19968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 19984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 20000(pA), rA0
            #if KB == 180
               prefB((pfB))
            #elif KB-2 == 180
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 180
               prefB(64(pfB))
            #elif KB-2 == 180
               prefA(64(pfA))
            #endif
            mulpd  20016(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 180
            MOVAPD 20032(pA), rA0
            MOVAPD 1312(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 20048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 20064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 20080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 20096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 20112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 20128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 20144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 20160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 20176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 20192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 20208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 20224(pA), rA0
            #if KB == 182
               prefB((pfB))
            #elif KB-2 == 182
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 182
               prefB(64(pfB))
            #elif KB-2 == 182
               prefA(64(pfA))
            #endif
            mulpd  20240(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 182
            MOVAPD 20256(pA), rA0
            MOVAPD 1328(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 20272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 20288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 20304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 20320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 20336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 20352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 20368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 20384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 20400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 20416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 20432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 20448(pA), rA0
            #if KB == 184
               prefB((pfB))
            #elif KB-2 == 184
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 184
               prefB(64(pfB))
            #elif KB-2 == 184
               prefA(64(pfA))
            #endif
            mulpd  20464(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 184
            MOVAPD 20480(pA), rA0
            MOVAPD 1344(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 20496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 20512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 20528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 20544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 20560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 20576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 20592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 20608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 20624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 20640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 20656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 20672(pA), rA0
            #if KB == 186
               prefB((pfB))
            #elif KB-2 == 186
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 186
               prefB(64(pfB))
            #elif KB-2 == 186
               prefA(64(pfA))
            #endif
            mulpd  20688(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 186
            MOVAPD 20704(pA), rA0
            MOVAPD 1360(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 20720(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 20736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 20752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 20768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 20784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 20800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 20816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 20832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 20848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 20864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 20880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 20896(pA), rA0
            #if KB == 188
               prefB((pfB))
            #elif KB-2 == 188
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 188
               prefB(64(pfB))
            #elif KB-2 == 188
               prefA(64(pfA))
            #endif
            mulpd  20912(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 188
            MOVAPD 20928(pA), rA0
            MOVAPD 1376(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 20944(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 20960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 20976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 20992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 21008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 21024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 21040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 21056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 21072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 21088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 21104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 21120(pA), rA0
            #if KB == 190
               prefB((pfB))
            #elif KB-2 == 190
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 190
               prefB(64(pfB))
            #elif KB-2 == 190
               prefA(64(pfA))
            #endif
            mulpd  21136(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 190
            MOVAPD 21152(pA), rA0
            MOVAPD 1392(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 21168(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 21184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 21200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 21216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 21232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 21248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 21264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 21280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 21296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 21312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 21328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 21344(pA), rA0
            #if KB == 192
               prefB((pfB))
            #elif KB-2 == 192
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 192
               prefB(64(pfB))
            #elif KB-2 == 192
               prefA(64(pfA))
            #endif
            mulpd  21360(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 192
            MOVAPD 21376(pA), rA0
            MOVAPD 1408(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 21392(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 21408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 21424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 21440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 21456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 21472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 21488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 21504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 21520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 21536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 21552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 21568(pA), rA0
            #if KB == 194
               prefB((pfB))
            #elif KB-2 == 194
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 194
               prefB(64(pfB))
            #elif KB-2 == 194
               prefA(64(pfA))
            #endif
            mulpd  21584(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 194
            MOVAPD 21600(pA), rA0
            MOVAPD 1424(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 21616(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 21632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 21648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 21664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 21680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 21696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 21712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 21728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 21744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 21760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 21776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 21792(pA), rA0
            #if KB == 196
               prefB((pfB))
            #elif KB-2 == 196
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 196
               prefB(64(pfB))
            #elif KB-2 == 196
               prefA(64(pfA))
            #endif
            mulpd  21808(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 196
            MOVAPD 21824(pA), rA0
            MOVAPD 1440(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 21840(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 21856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 21872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 21888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 21904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 21920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 21936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 21952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 21968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 21984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 22000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 22016(pA), rA0
            #if KB == 198
               prefB((pfB))
            #elif KB-2 == 198
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 198
               prefB(64(pfB))
            #elif KB-2 == 198
               prefA(64(pfA))
            #endif
            mulpd  22032(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 198
            MOVAPD 22048(pA), rA0
            MOVAPD 1456(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 22064(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 22080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 22096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 22112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 22128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 22144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 22160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 22176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 22192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 22208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 22224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 22240(pA), rA0
            #if KB == 200
               prefB((pfB))
            #elif KB-2 == 200
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 200
               prefB(64(pfB))
            #elif KB-2 == 200
               prefA(64(pfA))
            #endif
            mulpd  22256(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 200
            MOVAPD 22272(pA), rA0
            MOVAPD 1472(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 22288(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 22304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 22320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 22336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 22352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 22368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 22384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 22400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 22416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 22432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 22448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 22464(pA), rA0
            #if KB == 202
               prefB((pfB))
            #elif KB-2 == 202
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 202
               prefB(64(pfB))
            #elif KB-2 == 202
               prefA(64(pfA))
            #endif
            mulpd  22480(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 202
            MOVAPD 22496(pA), rA0
            MOVAPD 1488(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 22512(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 22528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 22544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 22560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 22576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 22592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 22608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 22624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 22640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 22656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 22672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 22688(pA), rA0
            #if KB == 204
               prefB((pfB))
            #elif KB-2 == 204
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 204
               prefB(64(pfB))
            #elif KB-2 == 204
               prefA(64(pfA))
            #endif
            mulpd  22704(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 204
            MOVAPD 22720(pA), rA0
            MOVAPD 1504(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 22736(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 22752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 22768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 22784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 22800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 22816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 22832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 22848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 22864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 22880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 22896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 22912(pA), rA0
            #if KB == 206
               prefB((pfB))
            #elif KB-2 == 206
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 206
               prefB(64(pfB))
            #elif KB-2 == 206
               prefA(64(pfA))
            #endif
            mulpd  22928(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 206
            MOVAPD 22944(pA), rA0
            MOVAPD 1520(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 22960(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 22976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 22992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 23008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 23024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 23040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 23056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 23072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 23088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 23104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 23120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 23136(pA), rA0
            #if KB == 208
               prefB((pfB))
            #elif KB-2 == 208
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 208
               prefB(64(pfB))
            #elif KB-2 == 208
               prefA(64(pfA))
            #endif
            mulpd  23152(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 208
            MOVAPD 23168(pA), rA0
            MOVAPD 1536(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 23184(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 23200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 23216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 23232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 23248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 23264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 23280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 23296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 23312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 23328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 23344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 23360(pA), rA0
            #if KB == 210
               prefB((pfB))
            #elif KB-2 == 210
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 210
               prefB(64(pfB))
            #elif KB-2 == 210
               prefA(64(pfA))
            #endif
            mulpd  23376(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 210
            MOVAPD 23392(pA), rA0
            MOVAPD 1552(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 23408(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 23424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 23440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 23456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 23472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 23488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 23504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 23520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 23536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 23552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 23568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 23584(pA), rA0
            #if KB == 212
               prefB((pfB))
            #elif KB-2 == 212
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 212
               prefB(64(pfB))
            #elif KB-2 == 212
               prefA(64(pfA))
            #endif
            mulpd  23600(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 212
            MOVAPD 23616(pA), rA0
            MOVAPD 1568(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 23632(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 23648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 23664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 23680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 23696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 23712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 23728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 23744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 23760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 23776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 23792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 23808(pA), rA0
            #if KB == 214
               prefB((pfB))
            #elif KB-2 == 214
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 214
               prefB(64(pfB))
            #elif KB-2 == 214
               prefA(64(pfA))
            #endif
            mulpd  23824(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 214
            MOVAPD 23840(pA), rA0
            MOVAPD 1584(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 23856(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 23872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 23888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 23904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 23920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 23936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 23952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 23968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 23984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 24000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 24016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 24032(pA), rA0
            #if KB == 216
               prefB((pfB))
            #elif KB-2 == 216
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 216
               prefB(64(pfB))
            #elif KB-2 == 216
               prefA(64(pfA))
            #endif
            mulpd  24048(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 216
            MOVAPD 24064(pA), rA0
            MOVAPD 1600(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 24080(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 24096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 24112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 24128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 24144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 24160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 24176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 24192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 24208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 24224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 24240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 24256(pA), rA0
            #if KB == 218
               prefB((pfB))
            #elif KB-2 == 218
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 218
               prefB(64(pfB))
            #elif KB-2 == 218
               prefA(64(pfA))
            #endif
            mulpd  24272(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 218
            MOVAPD 24288(pA), rA0
            MOVAPD 1616(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 24304(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 24320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 24336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 24352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 24368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 24384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 24400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 24416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 24432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 24448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 24464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 24480(pA), rA0
            #if KB == 220
               prefB((pfB))
            #elif KB-2 == 220
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 220
               prefB(64(pfB))
            #elif KB-2 == 220
               prefA(64(pfA))
            #endif
            mulpd  24496(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 220
            MOVAPD 24512(pA), rA0
            MOVAPD 1632(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 24528(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 24544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 24560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 24576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 24592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 24608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 24624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 24640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 24656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 24672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 24688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 24704(pA), rA0
            #if KB == 222
               prefB((pfB))
            #elif KB-2 == 222
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 222
               prefB(64(pfB))
            #elif KB-2 == 222
               prefA(64(pfA))
            #endif
            mulpd  24720(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 222
            MOVAPD 24736(pA), rA0
            MOVAPD 1648(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 24752(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 24768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 24784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 24800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 24816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 24832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 24848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 24864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 24880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 24896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 24912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 24928(pA), rA0
            #if KB == 224
               prefB((pfB))
            #elif KB-2 == 224
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 224
               prefB(64(pfB))
            #elif KB-2 == 224
               prefA(64(pfA))
            #endif
            mulpd  24944(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 224
            MOVAPD 24960(pA), rA0
            MOVAPD 1664(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 24976(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 24992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 25008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 25024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 25040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 25056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 25072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 25088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 25104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 25120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 25136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 25152(pA), rA0
            #if KB == 226
               prefB((pfB))
            #elif KB-2 == 226
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 226
               prefB(64(pfB))
            #elif KB-2 == 226
               prefA(64(pfA))
            #endif
            mulpd  25168(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 226
            MOVAPD 25184(pA), rA0
            MOVAPD 1680(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 25200(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 25216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 25232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 25248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 25264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 25280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 25296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 25312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 25328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 25344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 25360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 25376(pA), rA0
            #if KB == 228
               prefB((pfB))
            #elif KB-2 == 228
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 228
               prefB(64(pfB))
            #elif KB-2 == 228
               prefA(64(pfA))
            #endif
            mulpd  25392(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 228
            MOVAPD 25408(pA), rA0
            MOVAPD 1696(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 25424(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 25440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 25456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 25472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 25488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 25504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 25520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 25536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 25552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 25568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 25584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 25600(pA), rA0
            #if KB == 230
               prefB((pfB))
            #elif KB-2 == 230
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 230
               prefB(64(pfB))
            #elif KB-2 == 230
               prefA(64(pfA))
            #endif
            mulpd  25616(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 230
            MOVAPD 25632(pA), rA0
            MOVAPD 1712(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 25648(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 25664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 25680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 25696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 25712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 25728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 25744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 25760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 25776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 25792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 25808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 25824(pA), rA0
            #if KB == 232
               prefB((pfB))
            #elif KB-2 == 232
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 232
               prefB(64(pfB))
            #elif KB-2 == 232
               prefA(64(pfA))
            #endif
            mulpd  25840(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 232
            MOVAPD 25856(pA), rA0
            MOVAPD 1728(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 25872(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 25888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 25904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 25920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 25936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 25952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 25968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 25984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 26000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 26016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 26032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 26048(pA), rA0
            #if KB == 234
               prefB((pfB))
            #elif KB-2 == 234
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 234
               prefB(64(pfB))
            #elif KB-2 == 234
               prefA(64(pfA))
            #endif
            mulpd  26064(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 234
            MOVAPD 26080(pA), rA0
            MOVAPD 1744(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 26096(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 26112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 26128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 26144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 26160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 26176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 26192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 26208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 26224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 26240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 26256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 26272(pA), rA0
            #if KB == 236
               prefB((pfB))
            #elif KB-2 == 236
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 236
               prefB(64(pfB))
            #elif KB-2 == 236
               prefA(64(pfA))
            #endif
            mulpd  26288(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 236
            MOVAPD 26304(pA), rA0
            MOVAPD 1760(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 26320(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 26336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 26352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 26368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 26384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 26400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 26416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 26432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 26448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 26464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 26480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 26496(pA), rA0
            #if KB == 238
               prefB((pfB))
            #elif KB-2 == 238
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 238
               prefB(64(pfB))
            #elif KB-2 == 238
               prefA(64(pfA))
            #endif
            mulpd  26512(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 238
            MOVAPD 26528(pA), rA0
            MOVAPD 1776(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 26544(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 26560(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 26576(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 26592(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 26608(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 26624(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 26640(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 26656(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 26672(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 26688(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 26704(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 26720(pA), rA0
            #if KB == 240
               prefB((pfB))
            #elif KB-2 == 240
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 240
               prefB(64(pfB))
            #elif KB-2 == 240
               prefA(64(pfA))
            #endif
            mulpd  26736(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 240
            MOVAPD 26752(pA), rA0
            MOVAPD 1792(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 26768(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 26784(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 26800(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 26816(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 26832(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 26848(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 26864(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 26880(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 26896(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 26912(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 26928(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 26944(pA), rA0
            #if KB == 242
               prefB((pfB))
            #elif KB-2 == 242
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 242
               prefB(64(pfB))
            #elif KB-2 == 242
               prefA(64(pfA))
            #endif
            mulpd  26960(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 242
            MOVAPD 26976(pA), rA0
            MOVAPD 1808(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 26992(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 27008(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 27024(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 27040(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 27056(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 27072(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 27088(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 27104(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 27120(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 27136(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 27152(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 27168(pA), rA0
            #if KB == 244
               prefB((pfB))
            #elif KB-2 == 244
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 244
               prefB(64(pfB))
            #elif KB-2 == 244
               prefA(64(pfA))
            #endif
            mulpd  27184(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 244
            MOVAPD 27200(pA), rA0
            MOVAPD 1824(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 27216(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 27232(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 27248(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 27264(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 27280(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 27296(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 27312(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 27328(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 27344(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 27360(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 27376(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 27392(pA), rA0
            #if KB == 246
               prefB((pfB))
            #elif KB-2 == 246
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 246
               prefB(64(pfB))
            #elif KB-2 == 246
               prefA(64(pfA))
            #endif
            mulpd  27408(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 246
            MOVAPD 27424(pA), rA0
            MOVAPD 1840(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 27440(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 27456(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 27472(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 27488(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 27504(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 27520(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 27536(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 27552(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 27568(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 27584(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 27600(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 27616(pA), rA0
            #if KB == 248
               prefB((pfB))
            #elif KB-2 == 248
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 248
               prefB(64(pfB))
            #elif KB-2 == 248
               prefA(64(pfA))
            #endif
            mulpd  27632(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 248
            MOVAPD 27648(pA), rA0
            MOVAPD 1856(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 27664(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 27680(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 27696(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 27712(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 27728(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 27744(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 27760(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 27776(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 27792(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 27808(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 27824(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 27840(pA), rA0
            #if KB == 250
               prefB((pfB))
            #elif KB-2 == 250
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 250
               prefB(64(pfB))
            #elif KB-2 == 250
               prefA(64(pfA))
            #endif
            mulpd  27856(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 250
            MOVAPD 27872(pA), rA0
            MOVAPD 1872(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 27888(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 27904(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 27920(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 27936(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 27952(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 27968(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 27984(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 28000(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 28016(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 28032(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 28048(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 28064(pA), rA0
            #if KB == 252
               prefB((pfB))
            #elif KB-2 == 252
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 252
               prefB(64(pfB))
            #elif KB-2 == 252
               prefA(64(pfA))
            #endif
            mulpd  28080(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 252
            MOVAPD 28096(pA), rA0
            MOVAPD 1888(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 28112(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 28128(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 28144(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 28160(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 28176(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 28192(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 28208(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 28224(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 28240(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 28256(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 28272(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 28288(pA), rA0
            #if KB == 254
               prefB((pfB))
            #elif KB-2 == 254
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 254
               prefB(64(pfB))
            #elif KB-2 == 254
               prefA(64(pfA))
            #endif
            mulpd  28304(pA), rB0
            addpd rB0, rC13
         #endif
         #if KB > 254
            MOVAPD 28320(pA), rA0
            MOVAPD 1904(pB), rB0
            mulpd rB0, rA0
            addpd rA0, rC0
            MOVAPD 28336(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC1
            MOVAPD 28352(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC2
            MOVAPD 28368(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC3
            MOVAPD 28384(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC4
            MOVAPD 28400(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC5
            MOVAPD 28416(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC6
            MOVAPD 28432(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC7
            MOVAPD 28448(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC8
            MOVAPD 28464(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC9
            MOVAPD 28480(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC10
            MOVAPD 28496(pA), rA0
            mulpd rB0, rA0
            addpd rA0, rC11
            MOVAPD 28512(pA), rA0
            #if KB == 256
               prefB((pfB))
            #elif KB-2 == 256
               prefA((pfA))
            #endif
            mulpd rB0, rA0
            addpd rA0, rC12
            #if KB == 256
               prefB(64(pfB))
            #elif KB-2 == 256
               prefA(64(pfA))
            #endif
            mulpd  28528(pA), rB0
            addpd rB0, rC13
         #endif
KDONE:
/*
 *       Sum up rCx regs
 */
         haddpd rC1, rC0        /* rC0 = {c1ab,c0ab} */
         #ifdef VOP
            VOP (pC), rC0
         #endif
         MOVAPD rC0, (pC)
         haddpd rC3, rC2        /* rC2 = {c3ab,c2ab} */
         #ifdef VOP
            VOP 16(pC), rC2
         #endif
         MOVAPD rC2, 16(pC)
         haddpd rC5, rC4        /* rC4 = {c5ab,c4ab} */
         #ifdef VOP
            VOP 32(pC), rC4
         #endif
         MOVAPD rC4, 32(pC)
         haddpd rC7, rC6        /* rC6 = {c7ab,c6ab} */
         #ifdef VOP
            VOP 48(pC), rC6
         #endif
         MOVAPD rC6, 48(pC)
         haddpd rC9, rC8        /* rC8 = {c9ab,c8ab} */
         #ifdef VOP
            VOP 64(pC), rC8
         #endif
         MOVAPD rC8, 64(pC)
         haddpd rC11, rC10      /* rC10 = {c11ab,c10ab} */
         #ifdef VOP
            VOP 80(pC), rC10
         #endif
         MOVAPD rC10, 80(pC)
         haddpd rC13, rC12      /* rC12 = {c13ab,c12ab} */
         #ifdef VOP
            VOP 96(pC), rC12
         #endif
         MOVAPD rC12, 96(pC)

            add r112, pfB
            add r112, pfA
         add r112, pC
         add $KB*8, pB
         sub $1, NNU
      jnz MLOOP
      add incA, pA   /* pA += 14*K*sizeof */
      mov pB0, pB
      sub $1, NMU
      mov NNU0, NNU
   jnz MLOOP

DONE:
   movq (%rsp), %rbp
   movq 8(%rsp), %rbx
   movq 16(%rsp), %r12
   movq 24(%rsp), %r13
   movq 32(%rsp), %r14
   movq 40(%rsp), %r15
   add  $FSIZE, %rsp
   ret
