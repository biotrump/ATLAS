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
#define nmu     %rdi
#define nnu     %rsi
#define nnu0    %r10
#define pA      %rcx
#define pB      %rax
#define pC      %r9
#define pfA     %rbp
#define pB0     %r12
#define incPF   %rbx
#define pfB     %rdx
#define incAm   %r11

#define rA0     %xmm0
#define rB0     %xmm1
#define rB1     %xmm2
#define rB2     %xmm3
#define rB3     %xmm4
#define rC00    %xmm5
#define rC01    %xmm6
#define rC02    %xmm7
#define rC03    %xmm8
#define rC04    %xmm9
#define rC05    %xmm10
#define rC06    %xmm11
#define rC07    %xmm12
#define rC08    %xmm13
#ifndef prefA
   #if KB <= 72
      #define prefA prefetcht2
   #else
      #define prefA prefetcht2
   #endif
#endif
#ifndef prefB
   #if KB <= 72
      #define prefB prefetcht0
   #else
      #define prefB prefetcht2
   #endif
#endif
#ifndef prefC
   #ifdef ATL_3DNowXXX
      #define prefC prefetchw
   #else
      #define prefC prefetcht0
   #endif
#endif
#ifdef BETAN1
   #define BETCOP subps
#else
   #define BETCOP addps
#endif
#define FSIZE 6*8
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
     prefetcht0 (pA)
   sub $FSIZE, %rsp
     prefetcht0 (%r8)
   movq    %rbp, 0(%rsp)
     prefetcht0 64(pA)
   movq    %rbx, 8(%rsp)
     prefetcht0 128(pA)
   movq    %r12, 16(%rsp)
     prefetcht0 192(pA)
/*
 * Load paramaters
 */
   movq %r8, pB
     prefetcht0 256(pA)
   mov nnu, nnu0
     prefetcht0 320(pA)
   movq FSIZE+16(%rsp), pfB     /* pf = pBn */
     prefetcht0 384(pA)
   movq FSIZE+8(%rsp), pfA      /* pfB = pAn */
     prefetcht0 448(pA)
   mov $48*1*4, incPF   /* incPF = MU*NU*sizeof */
/*
 * Maximize are small-op size by adding 128 to ptrs
 */
   sub $-128, pA
     prefetcht0 512(pA)
   sub $-128, pB
     prefetcht0 576(pA)
   mov $KB*32*4, incAm           /* incAm = KB*MU*size */
   movq pB, pB0

   ALIGN16
   .local MNLOOP
   MNLOOP:
/*
 *       Peel 1st iteration for zeroing or rCxx & prefetch
 */
#if 0
         xorps rC00, rC00
         xorps rC01, rC01
         xorps rC02, rC02
         xorps rC03, rC03
         xorps rC04, rC04
         xorps rC05, rC05
         xorps rC06, rC06
         xorps rC07, rC07
#else
         movaps -128(pB), rB3    /* fmisc */
         pshufd $0x00, rB3, rB0  /* fadd */
         movaps -128(pA), rC00   /* fmul */

         mulps rB0, rC00         /* fmul */
         pshufd $0x55, rB3, rB1  /* fadd */
         movaps -112(pA), rC01   /* fmisc */

         mulps rB0, rC01         /* fmul */
         pshufd $0xAA, rB3, rB2  /* fadd */
         movaps -96(pA), rC02     /* fmisc */

         mulps rB0, rC02         /* fmul */
         shufps $0xFF, rB3, rB3  /* fadd */
         movaps -80(pA), rC03     /* fmisc */

         mulps rB0, rC03         /* fmul */
         movaps -64(pA), rC04     /* fadd */
         movaps -48(pA), rC05     /* fmisc */


         mulps rB0, rC04
         movaps -32(pA), rC06     /* fadd  */
         movaps -16(pA), rC07     /* fmisc */

         mulps rB0, rC05
         prefC (pC)
         prefC 64(pC)

         mulps rB0, rC06         /* fmul */
         prefB (pfB)
         prefA (pfA)
//         prefetcht2 64(pfB)

         mulps rB0, rC07         /* fmul */
//         prefetcht2 64(pfA)

         #if KB > 1
            movaps (pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 16(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 32(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 48(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 64(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 80(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 96(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 112(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 2
            movaps 128(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 144(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 160(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 176(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 192(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 208(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 224(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 240(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 3
            movaps 256(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 272(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 288(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 304(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 320(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 336(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 352(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            mulps 368(pA), rB3
            addps rB3, rC07
            #if KB > 4
               movaps -112(pB), rB3   /* fmisc */
            #endif
         #endif
#endif
/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         ALIGN16
         #if KB > 4
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 384(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 400(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 416(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 432(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 448(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 464(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 480(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 496(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 5
            movaps 512(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 528(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 544(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 560(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 576(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 592(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 608(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 624(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 6
            movaps 640(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 656(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 672(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 688(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 704(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 720(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 736(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 752(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 7
            movaps 768(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 784(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 800(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 816(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 832(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 848(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 864(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 880(pA), rA0
            mulps 880(pA), rB3
            addps rB3, rC07
            #if KB > 8
               movaps -96(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 8
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 896(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 912(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 928(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 944(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 960(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 976(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 992(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 1008(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 9
            movaps 1024(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 1040(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 1056(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 1072(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 1088(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 1104(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 1120(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 1136(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 10
            movaps 1152(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 1168(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 1184(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 1200(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 1216(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 1232(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 1248(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 1264(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 11
            movaps 1280(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 1296(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 1312(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 1328(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 1344(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 1360(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 1376(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 1392(pA), rA0
            mulps 1392(pA), rB3
            addps rB3, rC07
            #if KB > 12
               movaps -80(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 12
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 1408(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 1424(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 1440(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 1456(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 1472(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 1488(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 1504(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 1520(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 13
            movaps 1536(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 1552(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 1568(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 1584(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 1600(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 1616(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 1632(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 1648(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 14
            movaps 1664(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 1680(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 1696(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 1712(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 1728(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 1744(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 1760(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 1776(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 15
            movaps 1792(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 1808(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 1824(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 1840(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 1856(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 1872(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 1888(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 1904(pA), rA0
            mulps 1904(pA), rB3
            addps rB3, rC07
            #if KB > 16
               movaps -64(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 16
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 1920(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 1936(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 1952(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 1968(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 1984(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 2000(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 2016(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 2032(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 17
            movaps 2048(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 2064(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 2080(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 2096(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 2112(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 2128(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 2144(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 2160(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 18
            movaps 2176(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 2192(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 2208(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 2224(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 2240(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 2256(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 2272(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 2288(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 19
            movaps 2304(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 2320(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 2336(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 2352(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 2368(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 2384(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 2400(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 2416(pA), rA0
            mulps 2416(pA), rB3
            addps rB3, rC07
            #if KB > 20
               movaps -48(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 20
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 2432(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 2448(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 2464(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 2480(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 2496(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 2512(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 2528(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 2544(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 21
            movaps 2560(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 2576(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 2592(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 2608(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 2624(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 2640(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 2656(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 2672(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 22
            movaps 2688(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 2704(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 2720(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 2736(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 2752(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 2768(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 2784(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 2800(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 23
            movaps 2816(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 2832(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 2848(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 2864(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 2880(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 2896(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 2912(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 2928(pA), rA0
            mulps 2928(pA), rB3
            addps rB3, rC07
            #if KB > 24
               movaps -32(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 24
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 2944(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 2960(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 2976(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 2992(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 3008(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 3024(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 3040(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 3056(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 25
            movaps 3072(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 3088(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 3104(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 3120(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 3136(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 3152(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 3168(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 3184(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 26
            movaps 3200(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 3216(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 3232(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 3248(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 3264(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 3280(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 3296(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 3312(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 27
            movaps 3328(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 3344(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 3360(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 3376(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 3392(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 3408(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 3424(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 3440(pA), rA0
            mulps 3440(pA), rB3
            addps rB3, rC07
            #if KB > 28
               movaps -16(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 28
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 3456(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 3472(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 3488(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 3504(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 3520(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 3536(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 3552(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 3568(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 29
            movaps 3584(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 3600(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 3616(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 3632(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 3648(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 3664(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 3680(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 3696(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 30
            movaps 3712(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 3728(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 3744(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 3760(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 3776(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 3792(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 3808(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 3824(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 31
            movaps 3840(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 3856(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 3872(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 3888(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 3904(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 3920(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 3936(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 3952(pA), rA0
            mulps 3952(pA), rB3
            addps rB3, rC07
            #if KB > 32
               movaps 0(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 32
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 3968(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 3984(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 4000(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 4016(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 4032(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 4048(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 4064(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 4080(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 33
            movaps 4096(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 4112(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 4128(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 4144(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 4160(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 4176(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 4192(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 4208(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 34
            movaps 4224(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 4240(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 4256(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 4272(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 4288(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 4304(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 4320(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 4336(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 35
            movaps 4352(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 4368(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 4384(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 4400(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 4416(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 4432(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 4448(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 4464(pA), rA0
            mulps 4464(pA), rB3
            addps rB3, rC07
            #if KB > 36
               movaps 16(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 36
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 4480(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 4496(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 4512(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 4528(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 4544(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 4560(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 4576(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 4592(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 37
            movaps 4608(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 4624(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 4640(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 4656(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 4672(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 4688(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 4704(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 4720(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 38
            movaps 4736(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 4752(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 4768(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 4784(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 4800(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 4816(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 4832(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 4848(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 39
            movaps 4864(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 4880(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 4896(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 4912(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 4928(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 4944(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 4960(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 4976(pA), rA0
            mulps 4976(pA), rB3
            addps rB3, rC07
            #if KB > 40
               movaps 32(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 40
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 4992(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 5008(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 5024(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 5040(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 5056(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 5072(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 5088(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 5104(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 41
            movaps 5120(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 5136(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 5152(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 5168(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 5184(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 5200(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 5216(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 5232(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 42
            movaps 5248(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 5264(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 5280(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 5296(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 5312(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 5328(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 5344(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 5360(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 43
            movaps 5376(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 5392(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 5408(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 5424(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 5440(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 5456(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 5472(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 5488(pA), rA0
            mulps 5488(pA), rB3
            addps rB3, rC07
            #if KB > 44
               movaps 48(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 44
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 5504(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 5520(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 5536(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 5552(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 5568(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 5584(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 5600(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 5616(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 45
            movaps 5632(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 5648(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 5664(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 5680(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 5696(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 5712(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 5728(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 5744(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 46
            movaps 5760(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 5776(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 5792(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 5808(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 5824(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 5840(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 5856(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 5872(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 47
            movaps 5888(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 5904(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 5920(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 5936(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 5952(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 5968(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 5984(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 6000(pA), rA0
            mulps 6000(pA), rB3
            addps rB3, rC07
            #if KB > 48
               movaps 64(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 48
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 6016(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 6032(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 6048(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 6064(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 6080(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 6096(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 6112(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 6128(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 49
            movaps 6144(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 6160(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 6176(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 6192(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 6208(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 6224(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 6240(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 6256(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 50
            movaps 6272(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 6288(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 6304(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 6320(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 6336(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 6352(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 6368(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 6384(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 51
            movaps 6400(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 6416(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 6432(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 6448(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 6464(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 6480(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 6496(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 6512(pA), rA0
            mulps 6512(pA), rB3
            addps rB3, rC07
            #if KB > 52
               movaps 80(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 52
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 6528(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 6544(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 6560(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 6576(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 6592(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 6608(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 6624(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 6640(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 53
            movaps 6656(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 6672(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 6688(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 6704(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 6720(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 6736(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 6752(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 6768(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 54
            movaps 6784(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 6800(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 6816(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 6832(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 6848(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 6864(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 6880(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 6896(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 55
            movaps 6912(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 6928(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 6944(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 6960(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 6976(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 6992(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 7008(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 7024(pA), rA0
            mulps 7024(pA), rB3
            addps rB3, rC07
            #if KB > 56
               movaps 96(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 56
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 7040(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 7056(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 7072(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 7088(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 7104(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 7120(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 7136(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 7152(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 57
            movaps 7168(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 7184(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 7200(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 7216(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 7232(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 7248(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 7264(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 7280(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 58
            movaps 7296(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 7312(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 7328(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 7344(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 7360(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 7376(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 7392(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 7408(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 59
            movaps 7424(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 7440(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 7456(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 7472(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 7488(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 7504(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 7520(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 7536(pA), rA0
            mulps 7536(pA), rB3
            addps rB3, rC07
            #if KB > 60
               movaps 112(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 60
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 7552(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 7568(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 7584(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 7600(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 7616(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 7632(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 7648(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 7664(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 61
            movaps 7680(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 7696(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 7712(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 7728(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 7744(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 7760(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 7776(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 7792(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 62
            movaps 7808(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 7824(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 7840(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 7856(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 7872(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 7888(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 7904(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 7920(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 63
            movaps 7936(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 7952(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 7968(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 7984(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 8000(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 8016(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 8032(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 8048(pA), rA0
            mulps 8048(pA), rB3
            addps rB3, rC07
            #if KB > 64
               movaps 128(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 64
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 8064(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 8080(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 8096(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 8112(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 8128(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 8144(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 8160(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 8176(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 65
            movaps 8192(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 8208(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 8224(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 8240(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 8256(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 8272(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 8288(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 8304(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 66
            movaps 8320(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 8336(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 8352(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 8368(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 8384(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 8400(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 8416(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 8432(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 67
            movaps 8448(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 8464(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 8480(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 8496(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 8512(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 8528(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 8544(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 8560(pA), rA0
            mulps 8560(pA), rB3
            addps rB3, rC07
            #if KB > 68
               movaps 144(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 68
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 8576(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 8592(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 8608(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 8624(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 8640(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 8656(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 8672(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 8688(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 69
            movaps 8704(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 8720(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 8736(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 8752(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 8768(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 8784(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 8800(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 8816(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 70
            movaps 8832(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 8848(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 8864(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 8880(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 8896(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 8912(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 8928(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 8944(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 71
            movaps 8960(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 8976(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 8992(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 9008(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 9024(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 9040(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 9056(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 9072(pA), rA0
            mulps 9072(pA), rB3
            addps rB3, rC07
            #if KB > 72
               movaps 160(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 72
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 9088(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 9104(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 9120(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 9136(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 9152(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 9168(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 9184(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 9200(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 73
            movaps 9216(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 9232(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 9248(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 9264(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 9280(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 9296(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 9312(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 9328(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 74
            movaps 9344(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 9360(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 9376(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 9392(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 9408(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 9424(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 9440(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 9456(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 75
            movaps 9472(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 9488(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 9504(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 9520(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 9536(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 9552(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 9568(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 9584(pA), rA0
            mulps 9584(pA), rB3
            addps rB3, rC07
            #if KB > 76
               movaps 176(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 76
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 9600(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 9616(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 9632(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 9648(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 9664(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 9680(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 9696(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 9712(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 77
            movaps 9728(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 9744(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 9760(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 9776(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 9792(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 9808(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 9824(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 9840(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 78
            movaps 9856(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 9872(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 9888(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 9904(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 9920(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 9936(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 9952(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 9968(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 79
            movaps 9984(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 10000(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 10016(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 10032(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 10048(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 10064(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 10080(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 10096(pA), rA0
            mulps 10096(pA), rB3
            addps rB3, rC07
            #if KB > 80
               movaps 192(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 80
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 10112(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 10128(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 10144(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 10160(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 10176(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 10192(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 10208(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 10224(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 81
            movaps 10240(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 10256(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 10272(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 10288(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 10304(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 10320(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 10336(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 10352(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 82
            movaps 10368(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 10384(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 10400(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 10416(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 10432(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 10448(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 10464(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 10480(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 83
            movaps 10496(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 10512(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 10528(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 10544(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 10560(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 10576(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 10592(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 10608(pA), rA0
            mulps 10608(pA), rB3
            addps rB3, rC07
            #if KB > 84
               movaps 208(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 84
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 10624(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 10640(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 10656(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 10672(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 10688(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 10704(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 10720(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 10736(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 85
            movaps 10752(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 10768(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 10784(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 10800(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 10816(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 10832(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 10848(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 10864(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 86
            movaps 10880(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 10896(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 10912(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 10928(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 10944(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 10960(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 10976(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 10992(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 87
            movaps 11008(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 11024(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 11040(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 11056(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 11072(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 11088(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 11104(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 11120(pA), rA0
            mulps 11120(pA), rB3
            addps rB3, rC07
            #if KB > 88
               movaps 224(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 88
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 11136(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 11152(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 11168(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 11184(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 11200(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 11216(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 11232(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 11248(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 89
            movaps 11264(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 11280(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 11296(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 11312(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 11328(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 11344(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 11360(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 11376(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 90
            movaps 11392(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 11408(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 11424(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 11440(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 11456(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 11472(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 11488(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 11504(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 91
            movaps 11520(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 11536(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 11552(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 11568(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 11584(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 11600(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 11616(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 11632(pA), rA0
            mulps 11632(pA), rB3
            addps rB3, rC07
            #if KB > 92
               movaps 240(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 92
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 11648(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 11664(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 11680(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 11696(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 11712(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 11728(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 11744(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 11760(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 93
            movaps 11776(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 11792(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 11808(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 11824(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 11840(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 11856(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 11872(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 11888(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 94
            movaps 11904(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 11920(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 11936(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 11952(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 11968(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 11984(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 12000(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 12016(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 95
            movaps 12032(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 12048(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 12064(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 12080(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 12096(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 12112(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 12128(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 12144(pA), rA0
            mulps 12144(pA), rB3
            addps rB3, rC07
            #if KB > 96
               movaps 256(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 96
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 12160(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 12176(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 12192(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 12208(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 12224(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 12240(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 12256(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 12272(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 97
            movaps 12288(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 12304(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 12320(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 12336(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 12352(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 12368(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 12384(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 12400(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 98
            movaps 12416(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 12432(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 12448(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 12464(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 12480(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 12496(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 12512(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 12528(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 99
            movaps 12544(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 12560(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 12576(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 12592(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 12608(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 12624(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 12640(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 12656(pA), rA0
            mulps 12656(pA), rB3
            addps rB3, rC07
            #if KB > 100
               movaps 272(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 100
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 12672(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 12688(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 12704(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 12720(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 12736(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 12752(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 12768(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 12784(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 101
            movaps 12800(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 12816(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 12832(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 12848(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 12864(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 12880(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 12896(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 12912(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 102
            movaps 12928(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 12944(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 12960(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 12976(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 12992(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 13008(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 13024(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 13040(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 103
            movaps 13056(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 13072(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 13088(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 13104(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 13120(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 13136(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 13152(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 13168(pA), rA0
            mulps 13168(pA), rB3
            addps rB3, rC07
            #if KB > 104
               movaps 288(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 104
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 13184(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 13200(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 13216(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 13232(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 13248(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 13264(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 13280(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 13296(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 105
            movaps 13312(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 13328(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 13344(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 13360(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 13376(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 13392(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 13408(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 13424(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 106
            movaps 13440(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 13456(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 13472(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 13488(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 13504(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 13520(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 13536(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 13552(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 107
            movaps 13568(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 13584(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 13600(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 13616(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 13632(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 13648(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 13664(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 13680(pA), rA0
            mulps 13680(pA), rB3
            addps rB3, rC07
            #if KB > 108
               movaps 304(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 108
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 13696(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 13712(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 13728(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 13744(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 13760(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 13776(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 13792(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 13808(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 109
            movaps 13824(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 13840(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 13856(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 13872(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 13888(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 13904(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 13920(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 13936(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 110
            movaps 13952(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 13968(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 13984(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 14000(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 14016(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 14032(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 14048(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 14064(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 111
            movaps 14080(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 14096(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 14112(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 14128(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 14144(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 14160(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 14176(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 14192(pA), rA0
            mulps 14192(pA), rB3
            addps rB3, rC07
            #if KB > 112
               movaps 320(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 112
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 14208(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 14224(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 14240(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 14256(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 14272(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 14288(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 14304(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 14320(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 113
            movaps 14336(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 14352(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 14368(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 14384(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 14400(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 14416(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 14432(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 14448(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 114
            movaps 14464(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 14480(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 14496(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 14512(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 14528(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 14544(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 14560(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 14576(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 115
            movaps 14592(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 14608(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 14624(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 14640(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 14656(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 14672(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 14688(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 14704(pA), rA0
            mulps 14704(pA), rB3
            addps rB3, rC07
            #if KB > 116
               movaps 336(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 116
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 14720(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 14736(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 14752(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 14768(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 14784(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 14800(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 14816(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 14832(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 117
            movaps 14848(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 14864(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 14880(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 14896(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 14912(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 14928(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 14944(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 14960(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 118
            movaps 14976(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 14992(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 15008(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 15024(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 15040(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 15056(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 15072(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 15088(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 119
            movaps 15104(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 15120(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 15136(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 15152(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 15168(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 15184(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 15200(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 15216(pA), rA0
            mulps 15216(pA), rB3
            addps rB3, rC07
            #if KB > 120
               movaps 352(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 120
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 15232(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 15248(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 15264(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 15280(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 15296(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 15312(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 15328(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 15344(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 121
            movaps 15360(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 15376(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 15392(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 15408(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 15424(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 15440(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 15456(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 15472(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 122
            movaps 15488(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 15504(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 15520(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 15536(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 15552(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 15568(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 15584(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 15600(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 123
            movaps 15616(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 15632(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 15648(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 15664(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 15680(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 15696(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 15712(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 15728(pA), rA0
            mulps 15728(pA), rB3
            addps rB3, rC07
            #if KB > 124
               movaps 368(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 124
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 15744(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 15760(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 15776(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 15792(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 15808(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 15824(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 15840(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 15856(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 125
            movaps 15872(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 15888(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 15904(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 15920(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 15936(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 15952(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 15968(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 15984(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 126
            movaps 16000(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 16016(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 16032(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 16048(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 16064(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 16080(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 16096(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 16112(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 127
            movaps 16128(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 16144(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 16160(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 16176(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 16192(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 16208(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 16224(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 16240(pA), rA0
            mulps 16240(pA), rB3
            addps rB3, rC07
            #if KB > 128
               movaps 384(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 128
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 16256(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 16272(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 16288(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 16304(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 16320(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 16336(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 16352(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 16368(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 129
            movaps 16384(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 16400(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 16416(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 16432(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 16448(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 16464(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 16480(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 16496(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 130
            movaps 16512(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 16528(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 16544(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 16560(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 16576(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 16592(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 16608(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 16624(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 131
            movaps 16640(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 16656(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 16672(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 16688(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 16704(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 16720(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 16736(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 16752(pA), rA0
            mulps 16752(pA), rB3
            addps rB3, rC07
            #if KB > 132
               movaps 400(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 132
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 16768(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 16784(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 16800(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 16816(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 16832(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 16848(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 16864(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 16880(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 133
            movaps 16896(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 16912(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 16928(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 16944(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 16960(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 16976(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 16992(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 17008(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 134
            movaps 17024(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 17040(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 17056(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 17072(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 17088(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 17104(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 17120(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 17136(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 135
            movaps 17152(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 17168(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 17184(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 17200(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 17216(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 17232(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 17248(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 17264(pA), rA0
            mulps 17264(pA), rB3
            addps rB3, rC07
            #if KB > 136
               movaps 416(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 136
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 17280(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 17296(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 17312(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 17328(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 17344(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 17360(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 17376(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 17392(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 137
            movaps 17408(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 17424(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 17440(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 17456(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 17472(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 17488(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 17504(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 17520(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 138
            movaps 17536(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 17552(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 17568(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 17584(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 17600(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 17616(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 17632(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 17648(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 139
            movaps 17664(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 17680(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 17696(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 17712(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 17728(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 17744(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 17760(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 17776(pA), rA0
            mulps 17776(pA), rB3
            addps rB3, rC07
            #if KB > 140
               movaps 432(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 140
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 17792(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 17808(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 17824(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 17840(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 17856(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 17872(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 17888(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 17904(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 141
            movaps 17920(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 17936(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 17952(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 17968(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 17984(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 18000(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 18016(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 18032(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 142
            movaps 18048(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 18064(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 18080(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 18096(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 18112(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 18128(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 18144(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 18160(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 143
            movaps 18176(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 18192(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 18208(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 18224(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 18240(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 18256(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 18272(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 18288(pA), rA0
            mulps 18288(pA), rB3
            addps rB3, rC07
            #if KB > 144
               movaps 448(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 144
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 18304(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 18320(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 18336(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 18352(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 18368(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 18384(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 18400(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 18416(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 145
            movaps 18432(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 18448(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 18464(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 18480(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 18496(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 18512(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 18528(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 18544(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 146
            movaps 18560(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 18576(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 18592(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 18608(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 18624(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 18640(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 18656(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 18672(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 147
            movaps 18688(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 18704(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 18720(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 18736(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 18752(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 18768(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 18784(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 18800(pA), rA0
            mulps 18800(pA), rB3
            addps rB3, rC07
            #if KB > 148
               movaps 464(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 148
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 18816(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 18832(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 18848(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 18864(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 18880(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 18896(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 18912(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 18928(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 149
            movaps 18944(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 18960(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 18976(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 18992(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 19008(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 19024(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 19040(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 19056(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 150
            movaps 19072(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 19088(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 19104(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 19120(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 19136(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 19152(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 19168(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 19184(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 151
            movaps 19200(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 19216(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 19232(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 19248(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 19264(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 19280(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 19296(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 19312(pA), rA0
            mulps 19312(pA), rB3
            addps rB3, rC07
            #if KB > 152
               movaps 480(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 152
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 19328(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 19344(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 19360(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 19376(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 19392(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 19408(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 19424(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 19440(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 153
            movaps 19456(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 19472(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 19488(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 19504(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 19520(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 19536(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 19552(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 19568(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 154
            movaps 19584(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 19600(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 19616(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 19632(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 19648(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 19664(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 19680(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 19696(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 155
            movaps 19712(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 19728(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 19744(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 19760(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 19776(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 19792(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 19808(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 19824(pA), rA0
            mulps 19824(pA), rB3
            addps rB3, rC07
            #if KB > 156
               movaps 496(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 156
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 19840(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 19856(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 19872(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 19888(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 19904(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 19920(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 19936(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 19952(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 157
            movaps 19968(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 19984(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 20000(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 20016(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 20032(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 20048(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 20064(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 20080(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 158
            movaps 20096(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 20112(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 20128(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 20144(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 20160(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 20176(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 20192(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 20208(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 159
            movaps 20224(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 20240(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 20256(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 20272(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 20288(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 20304(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 20320(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 20336(pA), rA0
            mulps 20336(pA), rB3
            addps rB3, rC07
            #if KB > 160
               movaps 512(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 160
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 20352(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 20368(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 20384(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 20400(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 20416(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 20432(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 20448(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 20464(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 161
            movaps 20480(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 20496(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 20512(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 20528(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 20544(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 20560(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 20576(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 20592(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 162
            movaps 20608(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 20624(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 20640(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 20656(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 20672(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 20688(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 20704(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 20720(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 163
            movaps 20736(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 20752(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 20768(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 20784(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 20800(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 20816(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 20832(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 20848(pA), rA0
            mulps 20848(pA), rB3
            addps rB3, rC07
            #if KB > 164
               movaps 528(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 164
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 20864(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 20880(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 20896(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 20912(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 20928(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 20944(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 20960(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 20976(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 165
            movaps 20992(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 21008(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 21024(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 21040(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 21056(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 21072(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 21088(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 21104(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 166
            movaps 21120(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 21136(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 21152(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 21168(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 21184(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 21200(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 21216(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 21232(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 167
            movaps 21248(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 21264(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 21280(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 21296(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 21312(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 21328(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 21344(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 21360(pA), rA0
            mulps 21360(pA), rB3
            addps rB3, rC07
            #if KB > 168
               movaps 544(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 168
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 21376(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 21392(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 21408(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 21424(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 21440(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 21456(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 21472(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 21488(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 169
            movaps 21504(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 21520(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 21536(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 21552(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 21568(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 21584(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 21600(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 21616(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 170
            movaps 21632(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 21648(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 21664(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 21680(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 21696(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 21712(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 21728(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 21744(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 171
            movaps 21760(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 21776(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 21792(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 21808(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 21824(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 21840(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 21856(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 21872(pA), rA0
            mulps 21872(pA), rB3
            addps rB3, rC07
            #if KB > 172
               movaps 560(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 172
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 21888(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 21904(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 21920(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 21936(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 21952(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 21968(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 21984(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 22000(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 173
            movaps 22016(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 22032(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 22048(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 22064(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 22080(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 22096(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 22112(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 22128(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 174
            movaps 22144(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 22160(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 22176(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 22192(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 22208(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 22224(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 22240(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 22256(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 175
            movaps 22272(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 22288(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 22304(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 22320(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 22336(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 22352(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 22368(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 22384(pA), rA0
            mulps 22384(pA), rB3
            addps rB3, rC07
            #if KB > 176
               movaps 576(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 176
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 22400(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 22416(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 22432(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 22448(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 22464(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 22480(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 22496(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 22512(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 177
            movaps 22528(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 22544(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 22560(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 22576(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 22592(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 22608(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 22624(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 22640(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 178
            movaps 22656(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 22672(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 22688(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 22704(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 22720(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 22736(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 22752(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 22768(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 179
            movaps 22784(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 22800(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 22816(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 22832(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 22848(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 22864(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 22880(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 22896(pA), rA0
            mulps 22896(pA), rB3
            addps rB3, rC07
            #if KB > 180
               movaps 592(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 180
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 22912(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 22928(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 22944(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 22960(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 22976(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 22992(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 23008(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 23024(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 181
            movaps 23040(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 23056(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 23072(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 23088(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 23104(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 23120(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 23136(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 23152(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 182
            movaps 23168(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 23184(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 23200(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 23216(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 23232(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 23248(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 23264(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 23280(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 183
            movaps 23296(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 23312(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 23328(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 23344(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 23360(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 23376(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 23392(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 23408(pA), rA0
            mulps 23408(pA), rB3
            addps rB3, rC07
            #if KB > 184
               movaps 608(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 184
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 23424(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 23440(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 23456(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 23472(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 23488(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 23504(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 23520(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 23536(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 185
            movaps 23552(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 23568(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 23584(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 23600(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 23616(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 23632(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 23648(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 23664(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 186
            movaps 23680(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 23696(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 23712(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 23728(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 23744(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 23760(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 23776(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 23792(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 187
            movaps 23808(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 23824(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 23840(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 23856(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 23872(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 23888(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 23904(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 23920(pA), rA0
            mulps 23920(pA), rB3
            addps rB3, rC07
            #if KB > 188
               movaps 624(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 188
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 23936(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 23952(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 23968(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 23984(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 24000(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 24016(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 24032(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 24048(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 189
            movaps 24064(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 24080(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 24096(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 24112(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 24128(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 24144(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 24160(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 24176(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 190
            movaps 24192(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 24208(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 24224(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 24240(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 24256(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 24272(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 24288(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 24304(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 191
            movaps 24320(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 24336(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 24352(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 24368(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 24384(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 24400(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 24416(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 24432(pA), rA0
            mulps 24432(pA), rB3
            addps rB3, rC07
            #if KB > 192
               movaps 640(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 192
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 24448(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 24464(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 24480(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 24496(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 24512(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 24528(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 24544(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 24560(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 193
            movaps 24576(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 24592(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 24608(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 24624(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 24640(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 24656(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 24672(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 24688(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 194
            movaps 24704(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 24720(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 24736(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 24752(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 24768(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 24784(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 24800(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 24816(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 195
            movaps 24832(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 24848(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 24864(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 24880(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 24896(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 24912(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 24928(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 24944(pA), rA0
            mulps 24944(pA), rB3
            addps rB3, rC07
            #if KB > 196
               movaps 656(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 196
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 24960(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 24976(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 24992(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 25008(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 25024(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 25040(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 25056(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 25072(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 197
            movaps 25088(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 25104(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 25120(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 25136(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 25152(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 25168(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 25184(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 25200(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 198
            movaps 25216(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 25232(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 25248(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 25264(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 25280(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 25296(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 25312(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 25328(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 199
            movaps 25344(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 25360(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 25376(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 25392(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 25408(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 25424(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 25440(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 25456(pA), rA0
            mulps 25456(pA), rB3
            addps rB3, rC07
            #if KB > 200
               movaps 672(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 200
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 25472(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 25488(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 25504(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 25520(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 25536(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 25552(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 25568(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 25584(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 201
            movaps 25600(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 25616(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 25632(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 25648(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 25664(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 25680(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 25696(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 25712(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 202
            movaps 25728(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 25744(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 25760(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 25776(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 25792(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 25808(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 25824(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 25840(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 203
            movaps 25856(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 25872(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 25888(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 25904(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 25920(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 25936(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 25952(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 25968(pA), rA0
            mulps 25968(pA), rB3
            addps rB3, rC07
            #if KB > 204
               movaps 688(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 204
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 25984(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 26000(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 26016(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 26032(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 26048(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 26064(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 26080(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 26096(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 205
            movaps 26112(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 26128(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 26144(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 26160(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 26176(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 26192(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 26208(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 26224(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 206
            movaps 26240(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 26256(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 26272(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 26288(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 26304(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 26320(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 26336(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 26352(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 207
            movaps 26368(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 26384(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 26400(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 26416(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 26432(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 26448(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 26464(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 26480(pA), rA0
            mulps 26480(pA), rB3
            addps rB3, rC07
            #if KB > 208
               movaps 704(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 208
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 26496(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 26512(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 26528(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 26544(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 26560(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 26576(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 26592(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 26608(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 209
            movaps 26624(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 26640(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 26656(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 26672(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 26688(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 26704(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 26720(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 26736(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 210
            movaps 26752(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 26768(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 26784(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 26800(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 26816(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 26832(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 26848(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 26864(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 211
            movaps 26880(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 26896(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 26912(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 26928(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 26944(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 26960(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 26976(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 26992(pA), rA0
            mulps 26992(pA), rB3
            addps rB3, rC07
            #if KB > 212
               movaps 720(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 212
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 27008(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 27024(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 27040(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 27056(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 27072(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 27088(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 27104(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 27120(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 213
            movaps 27136(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 27152(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 27168(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 27184(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 27200(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 27216(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 27232(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 27248(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 214
            movaps 27264(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 27280(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 27296(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 27312(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 27328(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 27344(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 27360(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 27376(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 215
            movaps 27392(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 27408(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 27424(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 27440(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 27456(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 27472(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 27488(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 27504(pA), rA0
            mulps 27504(pA), rB3
            addps rB3, rC07
            #if KB > 216
               movaps 736(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 216
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 27520(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 27536(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 27552(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 27568(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 27584(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 27600(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 27616(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 27632(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 217
            movaps 27648(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 27664(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 27680(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 27696(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 27712(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 27728(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 27744(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 27760(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 218
            movaps 27776(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 27792(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 27808(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 27824(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 27840(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 27856(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 27872(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 27888(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 219
            movaps 27904(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 27920(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 27936(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 27952(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 27968(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 27984(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 28000(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 28016(pA), rA0
            mulps 28016(pA), rB3
            addps rB3, rC07
            #if KB > 220
               movaps 752(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 220
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 28032(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 28048(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 28064(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 28080(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 28096(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 28112(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 28128(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 28144(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 221
            movaps 28160(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 28176(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 28192(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 28208(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 28224(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 28240(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 28256(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 28272(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 222
            movaps 28288(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 28304(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 28320(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 28336(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 28352(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 28368(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 28384(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 28400(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 223
            movaps 28416(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 28432(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 28448(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 28464(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 28480(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 28496(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 28512(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 28528(pA), rA0
            mulps 28528(pA), rB3
            addps rB3, rC07
            #if KB > 224
               movaps 768(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 224
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 28544(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 28560(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 28576(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 28592(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 28608(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 28624(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 28640(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 28656(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 225
            movaps 28672(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 28688(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 28704(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 28720(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 28736(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 28752(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 28768(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 28784(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 226
            movaps 28800(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 28816(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 28832(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 28848(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 28864(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 28880(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 28896(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 28912(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 227
            movaps 28928(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 28944(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 28960(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 28976(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 28992(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 29008(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 29024(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 29040(pA), rA0
            mulps 29040(pA), rB3
            addps rB3, rC07
            #if KB > 228
               movaps 784(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 228
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 29056(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 29072(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 29088(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 29104(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 29120(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 29136(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 29152(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 29168(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 229
            movaps 29184(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 29200(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 29216(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 29232(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 29248(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 29264(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 29280(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 29296(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 230
            movaps 29312(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 29328(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 29344(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 29360(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 29376(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 29392(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 29408(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 29424(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 231
            movaps 29440(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 29456(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 29472(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 29488(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 29504(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 29520(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 29536(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 29552(pA), rA0
            mulps 29552(pA), rB3
            addps rB3, rC07
            #if KB > 232
               movaps 800(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 232
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 29568(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 29584(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 29600(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 29616(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 29632(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 29648(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 29664(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 29680(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 233
            movaps 29696(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 29712(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 29728(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 29744(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 29760(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 29776(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 29792(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 29808(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 234
            movaps 29824(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 29840(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 29856(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 29872(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 29888(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 29904(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 29920(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 29936(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 235
            movaps 29952(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 29968(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 29984(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 30000(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 30016(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 30032(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 30048(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 30064(pA), rA0
            mulps 30064(pA), rB3
            addps rB3, rC07
            #if KB > 236
               movaps 816(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 236
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 30080(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 30096(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 30112(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 30128(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 30144(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 30160(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 30176(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 30192(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 237
            movaps 30208(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 30224(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 30240(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 30256(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 30272(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 30288(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 30304(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 30320(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 238
            movaps 30336(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 30352(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 30368(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 30384(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 30400(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 30416(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 30432(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 30448(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 239
            movaps 30464(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 30480(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 30496(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 30512(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 30528(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 30544(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 30560(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 30576(pA), rA0
            mulps 30576(pA), rB3
            addps rB3, rC07
            #if KB > 240
               movaps 832(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 240
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 30592(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 30608(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 30624(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 30640(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 30656(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 30672(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 30688(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 30704(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 241
            movaps 30720(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 30736(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 30752(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 30768(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 30784(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 30800(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 30816(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 30832(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 242
            movaps 30848(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 30864(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 30880(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 30896(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 30912(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 30928(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 30944(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 30960(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 243
            movaps 30976(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 30992(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 31008(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 31024(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 31040(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 31056(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 31072(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 31088(pA), rA0
            mulps 31088(pA), rB3
            addps rB3, rC07
            #if KB > 244
               movaps 848(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 244
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 31104(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 31120(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 31136(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 31152(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 31168(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 31184(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 31200(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 31216(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 245
            movaps 31232(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 31248(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 31264(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 31280(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 31296(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 31312(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 31328(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 31344(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 246
            movaps 31360(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 31376(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 31392(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 31408(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 31424(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 31440(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 31456(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 31472(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 247
            movaps 31488(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 31504(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 31520(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 31536(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 31552(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 31568(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 31584(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 31600(pA), rA0
            mulps 31600(pA), rB3
            addps rB3, rC07
            #if KB > 248
               movaps 864(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 248
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 31616(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 31632(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 31648(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 31664(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 31680(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 31696(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 31712(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 31728(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 249
            movaps 31744(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 31760(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 31776(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 31792(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 31808(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 31824(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 31840(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 31856(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 250
            movaps 31872(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 31888(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 31904(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 31920(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 31936(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 31952(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 31968(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 31984(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 251
            movaps 32000(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 32016(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 32032(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 32048(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 32064(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 32080(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 32096(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 32112(pA), rA0
            mulps 32112(pA), rB3
            addps rB3, rC07
            #if KB > 252
               movaps 880(pB), rB3   /* fmisc */
            #endif
         #endif
         #if KB > 252
            pshufd $0x00, rB3, rB0  /* fadd */
            pshufd $0x55, rB3, rB1  /* fmul */
            movaps 32128(pA), rA0   /* fpmisc */
            pshufd $0xAA, rB3, rB2  /* fadd */
            shufps $0xFF, rB3, rB3  /* fmul */
            nop
            mulps rB0, rA0
            addps rA0, rC00
            movaps 32144(pA), rA0
            mulps rB0, rA0
            addps rA0, rC01
            movaps 32160(pA), rA0
            mulps rB0, rA0
            addps rA0, rC02
            movaps 32176(pA), rA0
            mulps rB0, rA0
            addps rA0, rC03
            movaps 32192(pA), rA0
            mulps rB0, rA0
            addps rA0, rC04
            movaps 32208(pA), rA0
            mulps rB0, rA0
            addps rA0, rC05
            movaps 32224(pA), rA0
            mulps rB0, rA0
            addps rA0, rC06
            mulps 32240(pA), rB0
            addps rB0, rC07
         #endif
         #if KB > 253
            movaps 32256(pA), rA0
            mulps rB1, rA0
            addps rA0, rC00
            movaps 32272(pA), rA0
            mulps rB1, rA0
            addps rA0, rC01
            movaps 32288(pA), rA0
            mulps rB1, rA0
            addps rA0, rC02
            movaps 32304(pA), rA0
            mulps rB1, rA0
            addps rA0, rC03
            movaps 32320(pA), rA0
            mulps rB1, rA0
            addps rA0, rC04
            movaps 32336(pA), rA0
            mulps rB1, rA0
            addps rA0, rC05
            movaps 32352(pA), rA0
            mulps rB1, rA0
            addps rA0, rC06
            mulps 32368(pA), rB1
            addps rB1, rC07
         #endif
         #if KB > 254
            movaps 32384(pA), rA0
            mulps rB2, rA0
            addps rA0, rC00
            movaps 32400(pA), rA0
            mulps rB2, rA0
            addps rA0, rC01
            movaps 32416(pA), rA0
            mulps rB2, rA0
            addps rA0, rC02
            movaps 32432(pA), rA0
            mulps rB2, rA0
            addps rA0, rC03
            movaps 32448(pA), rA0
            mulps rB2, rA0
            addps rA0, rC04
            movaps 32464(pA), rA0
            mulps rB2, rA0
            addps rA0, rC05
            movaps 32480(pA), rA0
            mulps rB2, rA0
            addps rA0, rC06
            mulps 32496(pA), rB2
            addps rB2, rC07
         #endif
         #if KB > 255
            movaps 32512(pA), rA0
            mulps rB3, rA0
            addps rA0, rC00
            movaps 32528(pA), rA0
            mulps rB3, rA0
            addps rA0, rC01
            movaps 32544(pA), rA0
            mulps rB3, rA0
            addps rA0, rC02
            movaps 32560(pA), rA0
            mulps rB3, rA0
            addps rA0, rC03
            movaps 32576(pA), rA0
            mulps rB3, rA0
            addps rA0, rC04
            movaps 32592(pA), rA0
            mulps rB3, rA0
            addps rA0, rC05
            movaps 32608(pA), rA0
            mulps rB3, rA0
            addps rA0, rC06
            movaps 32624(pA), rA0
            mulps 32624(pA), rB3
            addps rB3, rC07
            #if KB > 256
               movaps 896(pB), rB3   /* fmisc */
            #endif
         #endif
         #if defined(BETA1) || defined(BETAN1)
            BETCOP (pC), rC00
            movaps rC00, (pC)
            BETCOP 16(pC), rC01
            movaps rC01, 16(pC)
            BETCOP 32(pC), rC02
            movaps rC02, 32(pC)
            BETCOP 48(pC), rC03
            movaps rC03, 48(pC)
            BETCOP 64(pC), rC04
            movaps rC04, 64(pC)
            BETCOP 80(pC), rC05
            movaps rC05, 80(pC)
            BETCOP 96(pC), rC06
            movaps rC06, 96(pC)
            BETCOP 112(pC), rC07
            movaps rC07, 112(pC)
         #else
            movaps rC00, (pC)
            movaps rC01, 16(pC)
            movaps rC02, 32(pC)
            movaps rC03, 48(pC)
            movaps rC04, 64(pC)
            movaps rC05, 80(pC)
            movaps rC06, 96(pC)
            movaps rC07, 112(pC)
         #endif
         add incPF, pfA
         add incPF, pfB
         add $128, pC
         add $KB*1*4, pB   /* pB += KB*NU*sizeof */
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
   add $FSIZE, %rsp
   ret
#if 0
.global findSize
findSize:
mov $SS1-SS0, %rax
ret
SS0:
SS1:
#endif
