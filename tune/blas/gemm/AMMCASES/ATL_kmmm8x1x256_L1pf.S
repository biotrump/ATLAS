/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2015, 2008 R. Clint Whaley
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
 * This routine hastily adapted to the 3rd gen Opteron from the Core2 kernel,
 * ATL_dmm8x1x120.c.  Main difference is the prefetch strategy, which fetches
 * to the L1 on the Opt3rdGen, and the handling of loads, where MOVUPD
 * is preferred over MOVLPD/MOVHPD pair.  We also turn on unaligned mem
 * MULPD/ADDPD.  This all makes the kernal about 8% faster, and allows us
 * to use a smaller block factor (for better application performance).
 * This kernel does not do well for complex.
 * Later adapted to access-major in 2015
 */
#include "atlas_asm.h"
#ifndef ATL_GAS_x8664
   #error "This kernel requires x86-64 assembly!"
#endif
#if !defined(MB)
   #define MB 0
#endif
#if !defined(NB)
   #define NB 0
#endif
#if !defined(KB)
   #define KB 0
#endif
#if KB == 0
   #error "KB must be compile time constant!"
#endif
/*
 *Register usage
 */
#define pA0     %rcx
#define NN0     %rbx
#define pB      %rdx
#define pB0     %rax
#define pC0     %rsi
#define pfB     %rbp
#define pfA     %rdi
#define MM      %r10
#define NN      %r11

#define rA0     %xmm0
#define rB0     %xmm1
#define rC00    %xmm2
#define rC01    %xmm3
#define rC02    %xmm4
#define rC03    %xmm5
#define rC04    %xmm6
#define rC05    %xmm7
#define rC06    %xmm8
#define rC07    %xmm9
#define rC0     %xmm10
#define rC2     %xmm11
#define rC4     %xmm12
#define rC6     %xmm13

#if defined(SREAL) || defined(SCPLX)
   #define movapd movaps
   #define addpd addps
   #define subpd subps
   #define mulpd mulps
   #define SZ 4
   #define SZSH 2
   #define IMUL 2
#else
   #define SZ 8
   #define SZSH 3
   #define IMUL 1
#endif
#ifdef BETAN1
   #define VOP subpd
#elif defined(BETA1)
   #define VOP addpd
#elif defined(VOP)
   #undef VOP
#endif
/*
 * Prefetch defines
 */
#ifdef ATL_3DNow
   #if defined(ATL_UAMMID) && ATL_UAMMID == 0
      #define pref2(mem) prefetchw mem
   #else
      #define pref2(mem) prefetch mem
   #endif
   #define prefB(mem) prefetch  mem
   #define prefC(mem) prefetchw mem
#else
   #define pref2(mem) prefetcht0 mem
   #define prefB(mem) prefetcht0 mem
   #define prefC(mem) prefetcht0 mem
#endif
#if KB == 8
   #define KBSH 3
#elif KB == 16
   #define KBSH 4
#elif KB == 32
   #define KBSH 5
#elif KB == 64
   #define KBSH 6
#elif KB == 128
   #define KBSH 7
#elif defined(KBSH)
   #undef KBSH
#endif

#define PFAINC -8*SZ
#define PFBINC -8*SZ

/*
                    rdi      rsi    rdx        rcx         r8        r9
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                  8(%rsp)    16(%rsp)     24(%rsp)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */
        .text
.global ATL_asmdecor(ATL_USERMM)
ALIGN64
ATL_asmdecor(ATL_USERMM):
/*
 *  Save callee-saved iregs
 */
   #define FSIZE 16
   sub $FSIZE, %rsp
   movq %rbp, (%rsp)
   movq %rbx, 8(%rsp)
/*
 * Setup input parameters
 */
   mov  %rdi, MM
   mov  %rsi, NN
   mov  %r8, pB0
   mov  %r9, pC0
   mov MM, pfB
   #ifdef KBSH
      shl $KBSH+3+SZSH, pfB
   #else
      imul $KB*8*SZ, pfB
   #endif
   mov NN, pfA
   #ifdef KBSH
      shl $KBSH+SZSH, pfA
   #else
      imul $KB*SZ, pfA
   #endif
   mov FSIZE+16(%rsp), pB
   mov FSIZE+24(%rsp), NN0
   cmp pB, pB0
   cmovE NN0, pB
   add pB, pfA

//   addq FSIZE+16(%rsp), pfA
   addq FSIZE+8(%rsp), pfB
/*
 * pA0 += 128, pB0 += 128
 */
   sub     $-128, pA0
   sub     $-128, pB0
   mov pB0, pB
   mov NN, NN0
   .local MNLOOP
   ALIGN32
   MNLOOP:
            /*KLOOP: */
            movapd -128(pB0), rB0
                                        pref2((pfA))
                                        add     $PFAINC, pfA
            movapd -128(pA0), rC00
            mulpd rB0,rC00
            movapd -112(pA0), rC01
            mulpd rB0,rC01
            movapd -96(pA0), rC02
            mulpd rB0,rC02
            movapd -80(pA0), rC03
            mulpd rB0,rC03
            movapd -64(pA0), rC04
            mulpd rB0,rC04
            movapd -48(pA0), rC05
            mulpd rB0,rC05
            movapd -32(pA0), rC06
            mulpd rB0,rC06
            movapd -16(pA0), rC07
            mulpd rB0,rC07

            #if KB > 2*IMUL
               movapd -112(pB0), rB0
               movapd 0(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 16(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 32(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 48(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 64(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 80(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 96(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 112(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 4*IMUL
               movapd -96(pB0), rB0
               movapd 128(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 144(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 160(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 176(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 192(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 208(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 224(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 240(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 6*IMUL
               movapd -80(pB0), rB0
               movapd 256(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 272(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 288(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 304(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 320(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 336(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 352(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 368(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 8*IMUL
               movapd -64(pB0), rB0
               movapd 384(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 400(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 416(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 432(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 448(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 464(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 480(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 496(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 10*IMUL
               movapd -48(pB0), rB0
               movapd 512(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 528(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 544(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 560(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 576(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 592(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 608(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 624(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 12*IMUL
               movapd -32(pB0), rB0
               movapd 640(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 656(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 672(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 688(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 704(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 720(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 736(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 752(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 14*IMUL
               movapd -16(pB0), rB0
               movapd 768(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 784(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 800(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 816(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 832(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 848(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 864(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 880(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 16*IMUL
               movapd 0(pB0), rB0
               movapd 896(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 912(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 928(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 944(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 960(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 976(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 992(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 1008(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 18*IMUL
               movapd 16(pB0), rB0
               movapd 1024(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 1040(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 1056(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 1072(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 1088(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 1104(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 1120(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 1136(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 20*IMUL
               movapd 32(pB0), rB0
               movapd 1152(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 1168(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 1184(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 1200(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 1216(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 1232(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 1248(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 1264(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 22*IMUL
               movapd 48(pB0), rB0
               movapd 1280(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 1296(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 1312(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 1328(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 1344(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 1360(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 1376(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 1392(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 24*IMUL
               movapd 64(pB0), rB0
               movapd 1408(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 1424(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 1440(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 1456(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 1472(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 1488(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 1504(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 1520(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 26*IMUL
               movapd 80(pB0), rB0
               movapd 1536(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 1552(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 1568(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 1584(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 1600(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 1616(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 1632(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 1648(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 28*IMUL
               movapd 96(pB0), rB0
               movapd 1664(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 1680(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 1696(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 1712(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 1728(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 1744(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 1760(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 1776(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 30*IMUL
               movapd 112(pB0), rB0
               movapd 1792(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 1808(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 1824(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 1840(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 1856(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 1872(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 1888(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 1904(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 32*IMUL
               movapd 128(pB0), rB0
               movapd 1920(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 1936(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 1952(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 1968(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 1984(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 2000(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 2016(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 2032(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 34*IMUL
               movapd 144(pB0), rB0
               movapd 2048(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 2064(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 2080(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 2096(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 2112(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 2128(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 2144(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 2160(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 36*IMUL
               movapd 160(pB0), rB0
               movapd 2176(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 2192(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 2208(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 2224(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 2240(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 2256(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 2272(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 2288(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 38*IMUL
               movapd 176(pB0), rB0
               movapd 2304(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 2320(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 2336(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 2352(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 2368(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 2384(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 2400(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 2416(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 40*IMUL
               movapd 192(pB0), rB0
               movapd 2432(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 2448(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 2464(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 2480(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 2496(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 2512(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 2528(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 2544(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 42*IMUL
               movapd 208(pB0), rB0
               movapd 2560(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 2576(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 2592(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 2608(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 2624(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 2640(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 2656(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 2672(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 44*IMUL
               movapd 224(pB0), rB0
               movapd 2688(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 2704(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 2720(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 2736(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 2752(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 2768(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 2784(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 2800(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 46*IMUL
               movapd 240(pB0), rB0
               movapd 2816(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 2832(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 2848(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 2864(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 2880(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 2896(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 2912(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 2928(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 48*IMUL
               movapd 256(pB0), rB0
               movapd 2944(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 2960(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 2976(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 2992(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 3008(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 3024(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 3040(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 3056(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 50*IMUL
               movapd 272(pB0), rB0
               movapd 3072(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 3088(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 3104(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 3120(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 3136(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 3152(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 3168(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 3184(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 52*IMUL
               movapd 288(pB0), rB0
               movapd 3200(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 3216(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 3232(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 3248(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 3264(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 3280(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 3296(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 3312(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 54*IMUL
               movapd 304(pB0), rB0
               movapd 3328(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 3344(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 3360(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 3376(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 3392(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 3408(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 3424(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 3440(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 56*IMUL
               movapd 320(pB0), rB0
               movapd 3456(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 3472(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 3488(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 3504(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 3520(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 3536(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 3552(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 3568(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 58*IMUL
               movapd 336(pB0), rB0
               movapd 3584(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 3600(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 3616(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 3632(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 3648(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 3664(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 3680(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 3696(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 60*IMUL
               movapd 352(pB0), rB0
               movapd 3712(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 3728(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 3744(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 3760(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 3776(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 3792(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 3808(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 3824(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 62*IMUL
               movapd 368(pB0), rB0
               movapd 3840(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 3856(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 3872(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 3888(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 3904(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 3920(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 3936(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 3952(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 64*IMUL
               movapd 384(pB0), rB0
               movapd 3968(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 3984(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 4000(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 4016(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 4032(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 4048(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 4064(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 4080(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 66*IMUL
               movapd 400(pB0), rB0
               movapd 4096(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 4112(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 4128(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 4144(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 4160(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 4176(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 4192(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 4208(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 68*IMUL
               movapd 416(pB0), rB0
               movapd 4224(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 4240(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 4256(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 4272(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 4288(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 4304(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 4320(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 4336(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 70*IMUL
               movapd 432(pB0), rB0
               movapd 4352(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 4368(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 4384(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 4400(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 4416(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 4432(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 4448(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 4464(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 72*IMUL
               movapd 448(pB0), rB0
               movapd 4480(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 4496(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 4512(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 4528(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 4544(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 4560(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 4576(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 4592(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 74*IMUL
               movapd 464(pB0), rB0
               movapd 4608(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 4624(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 4640(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 4656(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 4672(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 4688(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 4704(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 4720(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 76*IMUL
               movapd 480(pB0), rB0
               movapd 4736(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 4752(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 4768(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 4784(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 4800(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 4816(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 4832(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 4848(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 78*IMUL
               movapd 496(pB0), rB0
               movapd 4864(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 4880(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 4896(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 4912(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 4928(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 4944(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 4960(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 4976(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 80*IMUL
               movapd 512(pB0), rB0
               movapd 4992(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 5008(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 5024(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 5040(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 5056(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 5072(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 5088(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 5104(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 82*IMUL
               movapd 528(pB0), rB0
               movapd 5120(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 5136(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 5152(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 5168(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 5184(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 5200(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 5216(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 5232(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 84*IMUL
               movapd 544(pB0), rB0
               movapd 5248(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 5264(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 5280(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 5296(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 5312(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 5328(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 5344(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 5360(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 86*IMUL
               movapd 560(pB0), rB0
               movapd 5376(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 5392(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 5408(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 5424(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 5440(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 5456(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 5472(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 5488(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 88*IMUL
               movapd 576(pB0), rB0
               movapd 5504(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 5520(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 5536(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 5552(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 5568(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 5584(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 5600(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 5616(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 90*IMUL
               movapd 592(pB0), rB0
               movapd 5632(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 5648(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 5664(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 5680(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 5696(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 5712(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 5728(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 5744(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 92*IMUL
               movapd 608(pB0), rB0
               movapd 5760(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 5776(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 5792(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 5808(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 5824(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 5840(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 5856(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 5872(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 94*IMUL
               movapd 624(pB0), rB0
               movapd 5888(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 5904(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 5920(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 5936(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 5952(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 5968(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 5984(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 6000(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 96*IMUL
               movapd 640(pB0), rB0
               movapd 6016(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 6032(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 6048(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 6064(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 6080(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 6096(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 6112(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 6128(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 98*IMUL
               movapd 656(pB0), rB0
               movapd 6144(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 6160(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 6176(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 6192(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 6208(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 6224(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 6240(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 6256(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 100*IMUL
               movapd 672(pB0), rB0
               movapd 6272(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 6288(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 6304(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 6320(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 6336(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 6352(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 6368(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 6384(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 102*IMUL
               movapd 688(pB0), rB0
               movapd 6400(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 6416(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 6432(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 6448(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 6464(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 6480(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 6496(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 6512(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 104*IMUL
               movapd 704(pB0), rB0
               movapd 6528(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 6544(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 6560(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 6576(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 6592(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 6608(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 6624(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 6640(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 106*IMUL
               movapd 720(pB0), rB0
               movapd 6656(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 6672(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 6688(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 6704(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 6720(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 6736(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 6752(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 6768(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 108*IMUL
               movapd 736(pB0), rB0
               movapd 6784(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 6800(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 6816(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 6832(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 6848(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 6864(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 6880(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 6896(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 110*IMUL
               movapd 752(pB0), rB0
               movapd 6912(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 6928(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 6944(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 6960(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 6976(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 6992(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 7008(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 7024(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 112*IMUL
               movapd 768(pB0), rB0
               movapd 7040(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 7056(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 7072(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 7088(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 7104(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 7120(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 7136(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 7152(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 114*IMUL
               movapd 784(pB0), rB0
               movapd 7168(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 7184(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 7200(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 7216(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 7232(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 7248(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 7264(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 7280(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 116*IMUL
               movapd 800(pB0), rB0
               movapd 7296(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 7312(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 7328(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 7344(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 7360(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 7376(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 7392(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 7408(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 118*IMUL
               movapd 816(pB0), rB0
               movapd 7424(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 7440(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 7456(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 7472(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 7488(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 7504(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 7520(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 7536(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 120*IMUL
               movapd 832(pB0), rB0
               movapd 7552(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 7568(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 7584(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 7600(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 7616(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 7632(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 7648(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 7664(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 122*IMUL
               movapd 848(pB0), rB0
               movapd 7680(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 7696(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 7712(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 7728(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 7744(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 7760(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 7776(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 7792(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 124*IMUL
               movapd 864(pB0), rB0
               movapd 7808(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 7824(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 7840(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 7856(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 7872(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 7888(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 7904(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 7920(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 126*IMUL
               movapd 880(pB0), rB0
               movapd 7936(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 7952(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 7968(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 7984(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 8000(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 8016(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 8032(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 8048(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 128*IMUL
               movapd 896(pB0), rB0
               movapd 8064(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 8080(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 8096(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 8112(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 8128(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 8144(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 8160(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 8176(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 130*IMUL
               movapd 912(pB0), rB0
               movapd 8192(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 8208(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 8224(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 8240(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 8256(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 8272(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 8288(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 8304(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 132*IMUL
               movapd 928(pB0), rB0
               movapd 8320(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 8336(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 8352(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 8368(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 8384(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 8400(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 8416(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 8432(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 134*IMUL
               movapd 944(pB0), rB0
               movapd 8448(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 8464(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 8480(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 8496(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 8512(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 8528(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 8544(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 8560(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 136*IMUL
               movapd 960(pB0), rB0
               movapd 8576(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 8592(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 8608(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 8624(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 8640(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 8656(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 8672(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 8688(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 138*IMUL
               movapd 976(pB0), rB0
               movapd 8704(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 8720(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 8736(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 8752(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 8768(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 8784(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 8800(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 8816(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 140*IMUL
               movapd 992(pB0), rB0
               movapd 8832(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 8848(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 8864(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 8880(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 8896(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 8912(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 8928(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 8944(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 142*IMUL
               movapd 1008(pB0), rB0
               movapd 8960(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 8976(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 8992(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 9008(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 9024(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 9040(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 9056(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 9072(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 144*IMUL
               movapd 1024(pB0), rB0
               movapd 9088(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 9104(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 9120(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 9136(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 9152(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 9168(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 9184(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 9200(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 146*IMUL
               movapd 1040(pB0), rB0
               movapd 9216(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 9232(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 9248(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 9264(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 9280(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 9296(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 9312(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 9328(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 148*IMUL
               movapd 1056(pB0), rB0
               movapd 9344(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 9360(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 9376(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 9392(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 9408(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 9424(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 9440(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 9456(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 150*IMUL
               movapd 1072(pB0), rB0
               movapd 9472(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 9488(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 9504(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 9520(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 9536(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 9552(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 9568(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 9584(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 152*IMUL
               movapd 1088(pB0), rB0
               movapd 9600(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 9616(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 9632(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 9648(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 9664(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 9680(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 9696(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 9712(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 154*IMUL
               movapd 1104(pB0), rB0
               movapd 9728(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 9744(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 9760(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 9776(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 9792(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 9808(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 9824(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 9840(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 156*IMUL
               movapd 1120(pB0), rB0
               movapd 9856(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 9872(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 9888(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 9904(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 9920(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 9936(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 9952(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 9968(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 158*IMUL
               movapd 1136(pB0), rB0
               movapd 9984(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 10000(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 10016(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 10032(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 10048(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 10064(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 10080(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 10096(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 160*IMUL
               movapd 1152(pB0), rB0
               movapd 10112(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 10128(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 10144(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 10160(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 10176(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 10192(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 10208(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 10224(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 162*IMUL
               movapd 1168(pB0), rB0
               movapd 10240(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 10256(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 10272(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 10288(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 10304(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 10320(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 10336(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 10352(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 164*IMUL
               movapd 1184(pB0), rB0
               movapd 10368(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 10384(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 10400(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 10416(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 10432(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 10448(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 10464(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 10480(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 166*IMUL
               movapd 1200(pB0), rB0
               movapd 10496(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 10512(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 10528(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 10544(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 10560(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 10576(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 10592(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 10608(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 168*IMUL
               movapd 1216(pB0), rB0
               movapd 10624(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 10640(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 10656(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 10672(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 10688(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 10704(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 10720(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 10736(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 170*IMUL
               movapd 1232(pB0), rB0
               movapd 10752(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 10768(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 10784(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 10800(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 10816(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 10832(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 10848(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 10864(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 172*IMUL
               movapd 1248(pB0), rB0
               movapd 10880(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 10896(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 10912(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 10928(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 10944(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 10960(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 10976(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 10992(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 174*IMUL
               movapd 1264(pB0), rB0
               movapd 11008(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 11024(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 11040(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 11056(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 11072(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 11088(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 11104(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 11120(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 176*IMUL
               movapd 1280(pB0), rB0
               movapd 11136(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 11152(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 11168(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 11184(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 11200(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 11216(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 11232(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 11248(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 178*IMUL
               movapd 1296(pB0), rB0
               movapd 11264(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 11280(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 11296(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 11312(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 11328(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 11344(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 11360(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 11376(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 180*IMUL
               movapd 1312(pB0), rB0
               movapd 11392(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 11408(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 11424(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 11440(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 11456(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 11472(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 11488(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 11504(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 182*IMUL
               movapd 1328(pB0), rB0
               movapd 11520(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 11536(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 11552(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 11568(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 11584(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 11600(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 11616(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 11632(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 184*IMUL
               movapd 1344(pB0), rB0
               movapd 11648(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 11664(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 11680(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 11696(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 11712(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 11728(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 11744(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 11760(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 186*IMUL
               movapd 1360(pB0), rB0
               movapd 11776(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 11792(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 11808(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 11824(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 11840(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 11856(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 11872(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 11888(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 188*IMUL
               movapd 1376(pB0), rB0
               movapd 11904(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 11920(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 11936(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 11952(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 11968(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 11984(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 12000(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 12016(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 190*IMUL
               movapd 1392(pB0), rB0
               movapd 12032(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 12048(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 12064(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 12080(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 12096(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 12112(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 12128(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 12144(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 192*IMUL
               movapd 1408(pB0), rB0
               movapd 12160(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 12176(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 12192(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 12208(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 12224(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 12240(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 12256(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 12272(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 194*IMUL
               movapd 1424(pB0), rB0
               movapd 12288(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 12304(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 12320(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 12336(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 12352(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 12368(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 12384(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 12400(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 196*IMUL
               movapd 1440(pB0), rB0
               movapd 12416(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 12432(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 12448(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 12464(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 12480(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 12496(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 12512(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 12528(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 198*IMUL
               movapd 1456(pB0), rB0
               movapd 12544(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 12560(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 12576(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 12592(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 12608(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 12624(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 12640(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 12656(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 200*IMUL
               movapd 1472(pB0), rB0
               movapd 12672(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 12688(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 12704(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 12720(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 12736(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 12752(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 12768(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 12784(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 202*IMUL
               movapd 1488(pB0), rB0
               movapd 12800(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 12816(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 12832(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 12848(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 12864(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 12880(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 12896(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 12912(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 204*IMUL
               movapd 1504(pB0), rB0
               movapd 12928(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 12944(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 12960(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 12976(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 12992(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 13008(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 13024(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 13040(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 206*IMUL
               movapd 1520(pB0), rB0
               movapd 13056(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 13072(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 13088(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 13104(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 13120(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 13136(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 13152(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 13168(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 208*IMUL
               movapd 1536(pB0), rB0
               movapd 13184(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 13200(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 13216(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 13232(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 13248(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 13264(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 13280(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 13296(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 210*IMUL
               movapd 1552(pB0), rB0
               movapd 13312(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 13328(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 13344(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 13360(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 13376(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 13392(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 13408(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 13424(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 212*IMUL
               movapd 1568(pB0), rB0
               movapd 13440(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 13456(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 13472(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 13488(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 13504(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 13520(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 13536(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 13552(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 214*IMUL
               movapd 1584(pB0), rB0
               movapd 13568(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 13584(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 13600(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 13616(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 13632(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 13648(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 13664(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 13680(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 216*IMUL
               movapd 1600(pB0), rB0
               movapd 13696(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 13712(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 13728(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 13744(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 13760(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 13776(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 13792(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 13808(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 218*IMUL
               movapd 1616(pB0), rB0
               movapd 13824(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 13840(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 13856(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 13872(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 13888(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 13904(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 13920(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 13936(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 220*IMUL
               movapd 1632(pB0), rB0
               movapd 13952(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 13968(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 13984(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 14000(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 14016(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 14032(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 14048(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 14064(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 222*IMUL
               movapd 1648(pB0), rB0
               movapd 14080(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 14096(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 14112(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 14128(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 14144(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 14160(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 14176(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 14192(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 224*IMUL
               movapd 1664(pB0), rB0
               movapd 14208(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 14224(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 14240(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 14256(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 14272(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 14288(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 14304(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 14320(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 226*IMUL
               movapd 1680(pB0), rB0
               movapd 14336(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 14352(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 14368(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 14384(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 14400(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 14416(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 14432(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 14448(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 228*IMUL
               movapd 1696(pB0), rB0
               movapd 14464(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 14480(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 14496(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 14512(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 14528(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 14544(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 14560(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 14576(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 230*IMUL
               movapd 1712(pB0), rB0
               movapd 14592(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 14608(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 14624(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 14640(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 14656(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 14672(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 14688(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 14704(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 232*IMUL
               movapd 1728(pB0), rB0
               movapd 14720(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 14736(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 14752(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 14768(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 14784(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 14800(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 14816(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 14832(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 234*IMUL
               movapd 1744(pB0), rB0
               movapd 14848(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 14864(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 14880(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 14896(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 14912(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 14928(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 14944(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 14960(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 236*IMUL
               movapd 1760(pB0), rB0
               movapd 14976(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 14992(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 15008(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 15024(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 15040(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 15056(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 15072(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 15088(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 238*IMUL
               movapd 1776(pB0), rB0
               movapd 15104(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 15120(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 15136(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 15152(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 15168(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 15184(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 15200(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 15216(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 240*IMUL
               movapd 1792(pB0), rB0
               movapd 15232(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 15248(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 15264(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 15280(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 15296(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 15312(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 15328(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 15344(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 242*IMUL
               movapd 1808(pB0), rB0
               movapd 15360(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 15376(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 15392(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 15408(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 15424(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 15440(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 15456(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 15472(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 244*IMUL
               movapd 1824(pB0), rB0
               movapd 15488(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 15504(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 15520(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 15536(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 15552(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 15568(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 15584(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 15600(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 246*IMUL
               movapd 1840(pB0), rB0
               movapd 15616(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 15632(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 15648(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 15664(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 15680(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 15696(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 15712(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 15728(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 248*IMUL
               movapd 1856(pB0), rB0
               movapd 15744(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 15760(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 15776(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 15792(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 15808(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 15824(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 15840(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 15856(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 250*IMUL
               movapd 1872(pB0), rB0
               movapd 15872(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 15888(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 15904(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 15920(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 15936(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 15952(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 15968(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 15984(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 252*IMUL
               movapd 1888(pB0), rB0
               movapd 16000(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 16016(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 16032(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 16048(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 16064(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 16080(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 16096(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 16112(pA0), rB0
               addpd rB0,rC07
            #endif
            #if KB > 254*IMUL
               movapd 1904(pB0), rB0
               movapd 16128(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC00
               movapd 16144(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC01
               movapd 16160(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC02
               movapd 16176(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC03
               movapd 16192(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC04
               movapd 16208(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC05
               movapd 16224(pA0), rA0
               mulpd rB0,rA0
               addpd rA0,rC06
               mulpd 16240(pA0), rB0
               addpd rB0,rC07
            #endif
            /* KLOOP end */
         /* rC00, rC01, rC02 .. rC07 */
         #if defined(SREAL) || defined(SCPLX)
            haddps rC01, rC00    /* rC00 = {c1cd,c1ab,c0cd,c0ab} */
            haddps rC03, rC02    /* rC20 = {c3cd,c3ab,c2cd,c2ab} */
            haddps rC02, rC00    /* rC00 = {c3, c2, c1, c0} */
            #ifdef VOP
               VOP (pC0), rC00
            #endif
            movaps rC00, (pC0)

            haddps rC05, rC04    /* rC04 = {c5cd,c5ab,c4cd,c4ab} */
            haddps rC07, rC06    /* rC26 = {c7cd,c7ab,c6cd,c6ab} */
            haddps rC06, rC04    /* rC04 = {c7, c6, c1, c0} */
            #ifdef VOP
               VOP 16(pC0), rC04
            #endif
            movaps rC04, 16(pC0)
         #elif !defined(ATL_SSE3) && defined(DREAL)
            movapd rC00, rC0     /* rC0  = {c0b, c0a} */
            movapd rC02, rC2     /* rC2  = {c2b, c2a} */
            unpcklpd rC01, rC00  /* rC00 = {c1a, c0a} */
            unpcklpd rC03, rC02  /* rC02 = {c3a, c2a} */
            unpckhpd rC01, rC0   /* rC0  = {c1b, c0b} */
            unpckhpd rC03, rC2   /* rC2  = {c3b, c2b} */
            addpd    rC0, rC00   /* rC00 = {c1ab, c0ab} */
            addpd    rC2, rC02   /* rC02 = {c3ab, c2ab} */

            movapd rC04, rC4     /* rC4  = {c4b, c4a} */
            movapd rC06, rC6     /* rC6  = {c6b, c6a} */
            unpcklpd rC05, rC04  /* rC04 = {c5a, c4a} */
            unpcklpd rC07, rC06  /* rC06 = {c7a, c6a} */
            unpckhpd rC05, rC4   /* rC4  = {c5b, c4b} */
            unpckhpd rC07, rC6   /* rC6  = {c7b, c6b} */
            addpd    rC4, rC04   /* rC04 = {c5ab, c4ab} */
            addpd    rC6, rC06   /* rC06 = {c7ab, c6ab} */
            #ifdef VOP
               VOP (pC0), rC00
               VOP 16(pC0), rC02
               VOP 32(pC0), rC04
               VOP 48(pC0), rC06
            #endif
            movapd rC00, (pC0)
            movapd rC02, 16(pC0)
            movapd rC04, 32(pC0)
            movapd rC06, 48(pC0)
         #else
/*
 *          pC[0-8] = rC[0-8]
 */
            haddpd  rC01,rC00
                  prefB((pfB))
                  add     $PFBINC, pfB
            #ifdef VOP
               VOP (pC0), rC00
            #endif
            movapd rC00, (pC0)

            haddpd  rC03,rC02
            #ifdef VOP
               VOP 16(pC0), rC02
            #endif
            movapd  rC02, 16(pC0)

            haddpd  rC05,rC04
            #ifdef VOP
               VOP 32(pC0), rC04
            #endif
            movapd  rC04, 32(pC0)

            haddpd  rC07,rC06
            #ifdef VOP
               VOP 48(pC0), rC06
            #endif
            movapd  rC06, 48(pC0)
         #endif

         add $KB*SZ, pB0
         add $8*SZ, pC0
         dec NN
      jnz MNLOOP
      add $8*KB*SZ, pA0
      mov pB, pB0
      movq NN0, NN
      dec MM
   jnz     MNLOOP

/* DONE: */
   movq (%rsp), %rbp
   movq 8(%rsp), %rbx
   add $FSIZE, %rsp
   ret
