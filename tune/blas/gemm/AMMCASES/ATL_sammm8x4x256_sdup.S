/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2014 R. Clint Whaley
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
#define movapd movaps
#define nmu     %rdi
#define nnu     %rsi
#define nnu0    %r10
#define pA      %rcx
#define pB      %rax
#define pC      %r9
#define pf      %rbp
#define pB0     %r12
#define incPF   %rbx
#define pfB     %rdx
#define incAm   %r11
#define pfC     %r13
#define r256    %r14

#define rm0     %xmm0
#define rA0     %xmm1
#define rA1     %xmm2
#define rB0     %xmm3
#define rB1     %xmm4
#define rB2     %xmm5
#define rB3     %xmm6
#define rC00    %xmm7
#define rC10    %xmm8
#define rC01    %xmm9
#define rC11    %xmm10
#define rC02    %xmm11
#define rC12    %xmm12
#define rC03    %xmm13
#define rC13    %xmm14
#ifndef pref
   #define pref prefetcht2
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
   #define BETCOP subps
#else
   #define BETCOP addps
#endif
#define FSIZE 4*8
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
/*
 * Load paramaters
 */
   movq %r8, pB
   mov nnu, nnu0
   movq FSIZE+16(%rsp), pf      /* pf = pBn */
   movq FSIZE+8(%rsp), pfB      /* pfB = pAn */
   movq FSIZE+24(%rsp), pfC    /* pfC = pCn */
   mov $8*8*4, incPF
/*
 * Extend range of small operands by starting at -128
 */
   sub $-128, pA
   sub $-128, pB
   mov $KB*8*4, incAm           /* incAm = KB*MU*size */
   movq pB, pB0

   ALIGN8
   .local MNLOOP
   MNLOOP:
/*
 *       Peel first iteration of K-loop to handle init of C to 0
 */
         movsldup -128(pB), rC02  /* port2, {b2, b2, b0, b0} */
         movddup rC02, rC00       /* port5, {b0, b0, b0, b0} */

         movaps -128(pA), rA0     /* port2, {a3, a2, a1, a0} */
         movaps rC00, rC10        /* port5 */
         mulps rA0, rC00          /* port0 */

         movaps -112(pA), rA1     /* port2 */
         mulps rA1, rC10          /* port0 */
         movhlps rC02, rC02       /* port5, {b2, b2, b2, b2} */

         movaps rC02, rC12        /* port5 */
         mulps rA0, rC02          /* port0 */
         movshdup -128(pB), rC03  /* port2, {b3, b3, b1, b1} */

         movddup rC03, rC01       /* port5, {b1, b1, b1, b1} */
         mulps rA1, rC12          /* port0 */
         #if KB > 1
            movsldup -112(pB), rB2/* port2, {b2, b2, b0, b0} */
         #else
            pref (pf)
         #endif

         movaps rC01, rC11        /* port5 */
         mulps rA0, rC01          /* port0 */
         #if KB > 1
            movshdup -112(pB), rB3/* port2, {b3, b3, b1, b1} */
         #else
            pref 64(pf)
         #endif

         movhlps rC03, rC03      /* port5, {b3, b3, b3, b3} */
         mulps rA1, rC11         /* port0 */

         movaps rC03, rC13       /* port5 */
         mulps rA0, rC03         /* port0 */

         #if KB > 1
            .byte 0x3e
            movddup rB2, rB0     /* port5 */
         #endif
         mulps rA1, rC13         /* port0 */
         #if KB == 1
            add incPF, pf        /* port[0,1,5] */
         #endif

/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 1
            movaps -96(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps -80(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 2
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 2
               movsldup -96(pB), rB2    /* port 2 */
            #elif KB == 2
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 2
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 2
               movshdup -96(pB), rB3    /* port 2 */
            #elif KB == 2
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 2
            movaps -64(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps -48(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 3
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 3
               movsldup -80(pB), rB2    /* port 2 */
            #elif KB == 3
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 3
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 3
               movshdup -80(pB), rB3    /* port 2 */
            #elif KB == 3
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 3
            movaps -32(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps -16(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 4
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 4
               movsldup -64(pB), rB2    /* port 2 */
            #elif KB == 4
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 4
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 4
               movshdup -64(pB), rB3    /* port 2 */
            #elif KB == 4
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 4
            movaps 0(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 16(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 5
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 5
               movsldup -48(pB), rB2    /* port 2 */
            #elif KB == 5
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 5
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 5
               movshdup -48(pB), rB3    /* port 2 */
            #elif KB == 5
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 5
            movaps 32(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 48(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 6
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 6
               movsldup -32(pB), rB2    /* port 2 */
            #elif KB == 6
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 6
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 6
               movshdup -32(pB), rB3    /* port 2 */
            #elif KB == 6
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 6
            movaps 64(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 80(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 7
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 7
               movsldup -16(pB), rB2    /* port 2 */
            #elif KB == 7
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 7
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 7
               movshdup -16(pB), rB3    /* port 2 */
            #elif KB == 7
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 7
            movaps 96(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 112(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 8
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 8
               movsldup 0(pB), rB2    /* port 2 */
            #elif KB == 8
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 8
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 8
               movshdup 0(pB), rB3    /* port 2 */
            #elif KB == 8
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 8
            movaps 128(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 144(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 9
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 9
               movsldup 16(pB), rB2    /* port 2 */
            #elif KB == 9
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 9
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 9
               movshdup 16(pB), rB3    /* port 2 */
            #elif KB == 9
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 9
            movaps 160(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 176(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 10
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 10
               movsldup 32(pB), rB2    /* port 2 */
            #elif KB == 10
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 10
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 10
               movshdup 32(pB), rB3    /* port 2 */
            #elif KB == 10
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 10
            movaps 192(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 208(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 11
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 11
               movsldup 48(pB), rB2    /* port 2 */
            #elif KB == 11
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 11
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 11
               movshdup 48(pB), rB3    /* port 2 */
            #elif KB == 11
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 11
            movaps 224(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 240(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 12
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 12
               movsldup 64(pB), rB2    /* port 2 */
            #elif KB == 12
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 12
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 12
               movshdup 64(pB), rB3    /* port 2 */
            #elif KB == 12
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 12
            movaps 256(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 272(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 13
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 13
               movsldup 80(pB), rB2    /* port 2 */
            #elif KB == 13
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 13
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 13
               movshdup 80(pB), rB3    /* port 2 */
            #elif KB == 13
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 13
            movaps 288(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 304(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 14
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 14
               movsldup 96(pB), rB2    /* port 2 */
            #elif KB == 14
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 14
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 14
               movshdup 96(pB), rB3    /* port 2 */
            #elif KB == 14
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 14
            movaps 320(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 336(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 15
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 15
               movsldup 112(pB), rB2    /* port 2 */
            #elif KB == 15
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 15
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 15
               movshdup 112(pB), rB3    /* port 2 */
            #elif KB == 15
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 15
            movaps 352(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 368(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 16
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 16
               movsldup 128(pB), rB2    /* port 2 */
            #elif KB == 16
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 16
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 16
               movshdup 128(pB), rB3    /* port 2 */
            #elif KB == 16
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 16
            movaps 384(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 400(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 17
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 17
               movsldup 144(pB), rB2    /* port 2 */
            #elif KB == 17
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 17
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 17
               movshdup 144(pB), rB3    /* port 2 */
            #elif KB == 17
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 17
            movaps 416(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 432(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 18
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 18
               movsldup 160(pB), rB2    /* port 2 */
            #elif KB == 18
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 18
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 18
               movshdup 160(pB), rB3    /* port 2 */
            #elif KB == 18
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 18
            movaps 448(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 464(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 19
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 19
               movsldup 176(pB), rB2    /* port 2 */
            #elif KB == 19
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 19
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 19
               movshdup 176(pB), rB3    /* port 2 */
            #elif KB == 19
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 19
            movaps 480(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 496(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 20
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 20
               movsldup 192(pB), rB2    /* port 2 */
            #elif KB == 20
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 20
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 20
               movshdup 192(pB), rB3    /* port 2 */
            #elif KB == 20
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 20
            movaps 512(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 528(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 21
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 21
               movsldup 208(pB), rB2    /* port 2 */
            #elif KB == 21
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 21
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 21
               movshdup 208(pB), rB3    /* port 2 */
            #elif KB == 21
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 21
            movaps 544(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 560(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 22
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 22
               movsldup 224(pB), rB2    /* port 2 */
            #elif KB == 22
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 22
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 22
               movshdup 224(pB), rB3    /* port 2 */
            #elif KB == 22
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 22
            movaps 576(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 592(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 23
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 23
               movsldup 240(pB), rB2    /* port 2 */
            #elif KB == 23
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 23
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 23
               movshdup 240(pB), rB3    /* port 2 */
            #elif KB == 23
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 23
            movaps 608(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 624(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 24
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 24
               movsldup 256(pB), rB2    /* port 2 */
            #elif KB == 24
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 24
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 24
               movshdup 256(pB), rB3    /* port 2 */
            #elif KB == 24
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 24
            movaps 640(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 656(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 25
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 25
               movsldup 272(pB), rB2    /* port 2 */
            #elif KB == 25
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 25
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 25
               movshdup 272(pB), rB3    /* port 2 */
            #elif KB == 25
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 25
            movaps 672(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 688(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 26
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 26
               movsldup 288(pB), rB2    /* port 2 */
            #elif KB == 26
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 26
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 26
               movshdup 288(pB), rB3    /* port 2 */
            #elif KB == 26
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 26
            movaps 704(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 720(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 27
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 27
               movsldup 304(pB), rB2    /* port 2 */
            #elif KB == 27
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 27
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 27
               movshdup 304(pB), rB3    /* port 2 */
            #elif KB == 27
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 27
            movaps 736(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 752(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 28
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 28
               movsldup 320(pB), rB2    /* port 2 */
            #elif KB == 28
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 28
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 28
               movshdup 320(pB), rB3    /* port 2 */
            #elif KB == 28
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 28
            movaps 768(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 784(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 29
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 29
               movsldup 336(pB), rB2    /* port 2 */
            #elif KB == 29
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 29
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 29
               movshdup 336(pB), rB3    /* port 2 */
            #elif KB == 29
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 29
            movaps 800(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 816(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 30
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 30
               movsldup 352(pB), rB2    /* port 2 */
            #elif KB == 30
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 30
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 30
               movshdup 352(pB), rB3    /* port 2 */
            #elif KB == 30
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 30
            movaps 832(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 848(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 31
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 31
               movsldup 368(pB), rB2    /* port 2 */
            #elif KB == 31
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 31
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 31
               movshdup 368(pB), rB3    /* port 2 */
            #elif KB == 31
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 31
            movaps 864(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 880(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 32
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 32
               movsldup 384(pB), rB2    /* port 2 */
            #elif KB == 32
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 32
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 32
               movshdup 384(pB), rB3    /* port 2 */
            #elif KB == 32
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 32
            movaps 896(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 912(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 33
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 33
               movsldup 400(pB), rB2    /* port 2 */
            #elif KB == 33
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 33
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 33
               movshdup 400(pB), rB3    /* port 2 */
            #elif KB == 33
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 33
            movaps 928(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 944(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 34
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 34
               movsldup 416(pB), rB2    /* port 2 */
            #elif KB == 34
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 34
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 34
               movshdup 416(pB), rB3    /* port 2 */
            #elif KB == 34
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 34
            movaps 960(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 976(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 35
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 35
               movsldup 432(pB), rB2    /* port 2 */
            #elif KB == 35
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 35
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 35
               movshdup 432(pB), rB3    /* port 2 */
            #elif KB == 35
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 35
            movaps 992(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1008(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 36
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 36
               movsldup 448(pB), rB2    /* port 2 */
            #elif KB == 36
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 36
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 36
               movshdup 448(pB), rB3    /* port 2 */
            #elif KB == 36
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 36
            movaps 1024(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1040(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 37
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 37
               movsldup 464(pB), rB2    /* port 2 */
            #elif KB == 37
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 37
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 37
               movshdup 464(pB), rB3    /* port 2 */
            #elif KB == 37
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 37
            movaps 1056(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1072(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 38
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 38
               movsldup 480(pB), rB2    /* port 2 */
            #elif KB == 38
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 38
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 38
               movshdup 480(pB), rB3    /* port 2 */
            #elif KB == 38
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 38
            movaps 1088(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1104(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 39
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 39
               movsldup 496(pB), rB2    /* port 2 */
            #elif KB == 39
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 39
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 39
               movshdup 496(pB), rB3    /* port 2 */
            #elif KB == 39
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 39
            movaps 1120(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1136(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 40
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 40
               movsldup 512(pB), rB2    /* port 2 */
            #elif KB == 40
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 40
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 40
               movshdup 512(pB), rB3    /* port 2 */
            #elif KB == 40
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 40
            movaps 1152(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1168(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 41
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 41
               movsldup 528(pB), rB2    /* port 2 */
            #elif KB == 41
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 41
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 41
               movshdup 528(pB), rB3    /* port 2 */
            #elif KB == 41
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 41
            movaps 1184(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1200(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 42
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 42
               movsldup 544(pB), rB2    /* port 2 */
            #elif KB == 42
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 42
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 42
               movshdup 544(pB), rB3    /* port 2 */
            #elif KB == 42
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 42
            movaps 1216(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1232(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 43
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 43
               movsldup 560(pB), rB2    /* port 2 */
            #elif KB == 43
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 43
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 43
               movshdup 560(pB), rB3    /* port 2 */
            #elif KB == 43
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 43
            movaps 1248(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1264(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 44
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 44
               movsldup 576(pB), rB2    /* port 2 */
            #elif KB == 44
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 44
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 44
               movshdup 576(pB), rB3    /* port 2 */
            #elif KB == 44
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 44
            movaps 1280(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1296(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 45
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 45
               movsldup 592(pB), rB2    /* port 2 */
            #elif KB == 45
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 45
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 45
               movshdup 592(pB), rB3    /* port 2 */
            #elif KB == 45
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 45
            movaps 1312(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1328(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 46
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 46
               movsldup 608(pB), rB2    /* port 2 */
            #elif KB == 46
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 46
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 46
               movshdup 608(pB), rB3    /* port 2 */
            #elif KB == 46
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 46
            movaps 1344(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1360(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 47
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 47
               movsldup 624(pB), rB2    /* port 2 */
            #elif KB == 47
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 47
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 47
               movshdup 624(pB), rB3    /* port 2 */
            #elif KB == 47
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 47
            movaps 1376(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1392(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 48
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 48
               movsldup 640(pB), rB2    /* port 2 */
            #elif KB == 48
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 48
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 48
               movshdup 640(pB), rB3    /* port 2 */
            #elif KB == 48
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 48
            movaps 1408(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1424(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 49
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 49
               movsldup 656(pB), rB2    /* port 2 */
            #elif KB == 49
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 49
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 49
               movshdup 656(pB), rB3    /* port 2 */
            #elif KB == 49
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 49
            movaps 1440(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1456(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 50
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 50
               movsldup 672(pB), rB2    /* port 2 */
            #elif KB == 50
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 50
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 50
               movshdup 672(pB), rB3    /* port 2 */
            #elif KB == 50
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 50
            movaps 1472(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1488(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 51
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 51
               movsldup 688(pB), rB2    /* port 2 */
            #elif KB == 51
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 51
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 51
               movshdup 688(pB), rB3    /* port 2 */
            #elif KB == 51
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 51
            movaps 1504(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1520(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 52
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 52
               movsldup 704(pB), rB2    /* port 2 */
            #elif KB == 52
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 52
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 52
               movshdup 704(pB), rB3    /* port 2 */
            #elif KB == 52
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 52
            movaps 1536(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1552(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 53
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 53
               movsldup 720(pB), rB2    /* port 2 */
            #elif KB == 53
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 53
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 53
               movshdup 720(pB), rB3    /* port 2 */
            #elif KB == 53
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 53
            movaps 1568(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1584(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 54
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 54
               movsldup 736(pB), rB2    /* port 2 */
            #elif KB == 54
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 54
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 54
               movshdup 736(pB), rB3    /* port 2 */
            #elif KB == 54
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 54
            movaps 1600(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1616(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 55
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 55
               movsldup 752(pB), rB2    /* port 2 */
            #elif KB == 55
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 55
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 55
               movshdup 752(pB), rB3    /* port 2 */
            #elif KB == 55
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 55
            movaps 1632(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1648(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 56
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 56
               movsldup 768(pB), rB2    /* port 2 */
            #elif KB == 56
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 56
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 56
               movshdup 768(pB), rB3    /* port 2 */
            #elif KB == 56
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 56
            movaps 1664(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1680(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 57
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 57
               movsldup 784(pB), rB2    /* port 2 */
            #elif KB == 57
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 57
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 57
               movshdup 784(pB), rB3    /* port 2 */
            #elif KB == 57
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 57
            movaps 1696(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1712(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 58
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 58
               movsldup 800(pB), rB2    /* port 2 */
            #elif KB == 58
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 58
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 58
               movshdup 800(pB), rB3    /* port 2 */
            #elif KB == 58
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 58
            movaps 1728(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1744(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 59
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 59
               movsldup 816(pB), rB2    /* port 2 */
            #elif KB == 59
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 59
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 59
               movshdup 816(pB), rB3    /* port 2 */
            #elif KB == 59
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 59
            movaps 1760(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1776(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 60
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 60
               movsldup 832(pB), rB2    /* port 2 */
            #elif KB == 60
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 60
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 60
               movshdup 832(pB), rB3    /* port 2 */
            #elif KB == 60
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 60
            movaps 1792(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1808(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 61
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 61
               movsldup 848(pB), rB2    /* port 2 */
            #elif KB == 61
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 61
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 61
               movshdup 848(pB), rB3    /* port 2 */
            #elif KB == 61
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 61
            movaps 1824(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1840(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 62
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 62
               movsldup 864(pB), rB2    /* port 2 */
            #elif KB == 62
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 62
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 62
               movshdup 864(pB), rB3    /* port 2 */
            #elif KB == 62
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 62
            movaps 1856(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1872(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 63
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 63
               movsldup 880(pB), rB2    /* port 2 */
            #elif KB == 63
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 63
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 63
               movshdup 880(pB), rB3    /* port 2 */
            #elif KB == 63
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 63
            movaps 1888(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1904(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 64
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 64
               movsldup 896(pB), rB2    /* port 2 */
            #elif KB == 64
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 64
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 64
               movshdup 896(pB), rB3    /* port 2 */
            #elif KB == 64
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 64
            movaps 1920(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1936(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 65
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 65
               movsldup 912(pB), rB2    /* port 2 */
            #elif KB == 65
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 65
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 65
               movshdup 912(pB), rB3    /* port 2 */
            #elif KB == 65
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 65
            movaps 1952(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1968(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 66
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 66
               movsldup 928(pB), rB2    /* port 2 */
            #elif KB == 66
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 66
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 66
               movshdup 928(pB), rB3    /* port 2 */
            #elif KB == 66
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 66
            movaps 1984(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2000(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 67
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 67
               movsldup 944(pB), rB2    /* port 2 */
            #elif KB == 67
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 67
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 67
               movshdup 944(pB), rB3    /* port 2 */
            #elif KB == 67
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 67
            movaps 2016(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2032(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 68
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 68
               movsldup 960(pB), rB2    /* port 2 */
            #elif KB == 68
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 68
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 68
               movshdup 960(pB), rB3    /* port 2 */
            #elif KB == 68
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 68
            movaps 2048(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2064(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 69
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 69
               movsldup 976(pB), rB2    /* port 2 */
            #elif KB == 69
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 69
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 69
               movshdup 976(pB), rB3    /* port 2 */
            #elif KB == 69
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 69
            movaps 2080(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2096(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 70
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 70
               movsldup 992(pB), rB2    /* port 2 */
            #elif KB == 70
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 70
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 70
               movshdup 992(pB), rB3    /* port 2 */
            #elif KB == 70
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 70
            movaps 2112(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2128(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 71
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 71
               movsldup 1008(pB), rB2    /* port 2 */
            #elif KB == 71
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 71
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 71
               movshdup 1008(pB), rB3    /* port 2 */
            #elif KB == 71
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 71
            movaps 2144(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2160(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 72
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 72
               movsldup 1024(pB), rB2    /* port 2 */
            #elif KB == 72
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 72
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 72
               movshdup 1024(pB), rB3    /* port 2 */
            #elif KB == 72
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 72
            movaps 2176(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2192(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 73
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 73
               movsldup 1040(pB), rB2    /* port 2 */
            #elif KB == 73
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 73
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 73
               movshdup 1040(pB), rB3    /* port 2 */
            #elif KB == 73
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 73
            movaps 2208(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2224(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 74
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 74
               movsldup 1056(pB), rB2    /* port 2 */
            #elif KB == 74
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 74
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 74
               movshdup 1056(pB), rB3    /* port 2 */
            #elif KB == 74
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 74
            movaps 2240(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2256(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 75
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 75
               movsldup 1072(pB), rB2    /* port 2 */
            #elif KB == 75
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 75
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 75
               movshdup 1072(pB), rB3    /* port 2 */
            #elif KB == 75
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 75
            movaps 2272(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2288(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 76
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 76
               movsldup 1088(pB), rB2    /* port 2 */
            #elif KB == 76
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 76
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 76
               movshdup 1088(pB), rB3    /* port 2 */
            #elif KB == 76
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 76
            movaps 2304(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2320(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 77
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 77
               movsldup 1104(pB), rB2    /* port 2 */
            #elif KB == 77
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 77
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 77
               movshdup 1104(pB), rB3    /* port 2 */
            #elif KB == 77
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 77
            movaps 2336(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2352(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 78
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 78
               movsldup 1120(pB), rB2    /* port 2 */
            #elif KB == 78
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 78
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 78
               movshdup 1120(pB), rB3    /* port 2 */
            #elif KB == 78
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 78
            movaps 2368(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2384(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 79
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 79
               movsldup 1136(pB), rB2    /* port 2 */
            #elif KB == 79
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 79
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 79
               movshdup 1136(pB), rB3    /* port 2 */
            #elif KB == 79
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 79
            movaps 2400(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2416(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 80
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 80
               movsldup 1152(pB), rB2    /* port 2 */
            #elif KB == 80
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 80
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 80
               movshdup 1152(pB), rB3    /* port 2 */
            #elif KB == 80
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 80
            movaps 2432(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2448(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 81
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 81
               movsldup 1168(pB), rB2    /* port 2 */
            #elif KB == 81
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 81
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 81
               movshdup 1168(pB), rB3    /* port 2 */
            #elif KB == 81
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 81
            movaps 2464(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2480(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 82
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 82
               movsldup 1184(pB), rB2    /* port 2 */
            #elif KB == 82
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 82
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 82
               movshdup 1184(pB), rB3    /* port 2 */
            #elif KB == 82
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 82
            movaps 2496(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2512(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 83
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 83
               movsldup 1200(pB), rB2    /* port 2 */
            #elif KB == 83
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 83
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 83
               movshdup 1200(pB), rB3    /* port 2 */
            #elif KB == 83
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 83
            movaps 2528(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2544(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 84
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 84
               movsldup 1216(pB), rB2    /* port 2 */
            #elif KB == 84
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 84
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 84
               movshdup 1216(pB), rB3    /* port 2 */
            #elif KB == 84
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 84
            movaps 2560(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2576(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 85
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 85
               movsldup 1232(pB), rB2    /* port 2 */
            #elif KB == 85
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 85
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 85
               movshdup 1232(pB), rB3    /* port 2 */
            #elif KB == 85
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 85
            movaps 2592(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2608(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 86
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 86
               movsldup 1248(pB), rB2    /* port 2 */
            #elif KB == 86
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 86
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 86
               movshdup 1248(pB), rB3    /* port 2 */
            #elif KB == 86
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 86
            movaps 2624(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2640(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 87
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 87
               movsldup 1264(pB), rB2    /* port 2 */
            #elif KB == 87
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 87
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 87
               movshdup 1264(pB), rB3    /* port 2 */
            #elif KB == 87
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 87
            movaps 2656(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2672(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 88
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 88
               movsldup 1280(pB), rB2    /* port 2 */
            #elif KB == 88
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 88
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 88
               movshdup 1280(pB), rB3    /* port 2 */
            #elif KB == 88
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 88
            movaps 2688(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2704(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 89
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 89
               movsldup 1296(pB), rB2    /* port 2 */
            #elif KB == 89
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 89
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 89
               movshdup 1296(pB), rB3    /* port 2 */
            #elif KB == 89
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 89
            movaps 2720(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2736(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 90
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 90
               movsldup 1312(pB), rB2    /* port 2 */
            #elif KB == 90
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 90
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 90
               movshdup 1312(pB), rB3    /* port 2 */
            #elif KB == 90
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 90
            movaps 2752(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2768(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 91
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 91
               movsldup 1328(pB), rB2    /* port 2 */
            #elif KB == 91
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 91
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 91
               movshdup 1328(pB), rB3    /* port 2 */
            #elif KB == 91
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 91
            movaps 2784(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2800(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 92
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 92
               movsldup 1344(pB), rB2    /* port 2 */
            #elif KB == 92
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 92
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 92
               movshdup 1344(pB), rB3    /* port 2 */
            #elif KB == 92
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 92
            movaps 2816(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2832(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 93
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 93
               movsldup 1360(pB), rB2    /* port 2 */
            #elif KB == 93
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 93
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 93
               movshdup 1360(pB), rB3    /* port 2 */
            #elif KB == 93
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 93
            movaps 2848(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2864(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 94
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 94
               movsldup 1376(pB), rB2    /* port 2 */
            #elif KB == 94
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 94
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 94
               movshdup 1376(pB), rB3    /* port 2 */
            #elif KB == 94
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 94
            movaps 2880(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2896(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 95
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 95
               movsldup 1392(pB), rB2    /* port 2 */
            #elif KB == 95
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 95
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 95
               movshdup 1392(pB), rB3    /* port 2 */
            #elif KB == 95
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 95
            movaps 2912(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2928(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 96
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 96
               movsldup 1408(pB), rB2    /* port 2 */
            #elif KB == 96
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 96
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 96
               movshdup 1408(pB), rB3    /* port 2 */
            #elif KB == 96
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 96
            movaps 2944(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2960(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 97
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 97
               movsldup 1424(pB), rB2    /* port 2 */
            #elif KB == 97
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 97
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 97
               movshdup 1424(pB), rB3    /* port 2 */
            #elif KB == 97
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 97
            movaps 2976(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2992(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 98
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 98
               movsldup 1440(pB), rB2    /* port 2 */
            #elif KB == 98
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 98
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 98
               movshdup 1440(pB), rB3    /* port 2 */
            #elif KB == 98
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 98
            movaps 3008(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3024(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 99
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 99
               movsldup 1456(pB), rB2    /* port 2 */
            #elif KB == 99
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 99
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 99
               movshdup 1456(pB), rB3    /* port 2 */
            #elif KB == 99
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 99
            movaps 3040(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3056(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 100
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 100
               movsldup 1472(pB), rB2    /* port 2 */
            #elif KB == 100
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 100
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 100
               movshdup 1472(pB), rB3    /* port 2 */
            #elif KB == 100
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 100
            movaps 3072(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3088(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 101
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 101
               movsldup 1488(pB), rB2    /* port 2 */
            #elif KB == 101
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 101
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 101
               movshdup 1488(pB), rB3    /* port 2 */
            #elif KB == 101
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 101
            movaps 3104(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3120(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 102
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 102
               movsldup 1504(pB), rB2    /* port 2 */
            #elif KB == 102
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 102
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 102
               movshdup 1504(pB), rB3    /* port 2 */
            #elif KB == 102
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 102
            movaps 3136(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3152(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 103
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 103
               movsldup 1520(pB), rB2    /* port 2 */
            #elif KB == 103
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 103
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 103
               movshdup 1520(pB), rB3    /* port 2 */
            #elif KB == 103
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 103
            movaps 3168(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3184(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 104
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 104
               movsldup 1536(pB), rB2    /* port 2 */
            #elif KB == 104
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 104
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 104
               movshdup 1536(pB), rB3    /* port 2 */
            #elif KB == 104
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 104
            movaps 3200(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3216(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 105
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 105
               movsldup 1552(pB), rB2    /* port 2 */
            #elif KB == 105
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 105
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 105
               movshdup 1552(pB), rB3    /* port 2 */
            #elif KB == 105
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 105
            movaps 3232(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3248(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 106
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 106
               movsldup 1568(pB), rB2    /* port 2 */
            #elif KB == 106
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 106
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 106
               movshdup 1568(pB), rB3    /* port 2 */
            #elif KB == 106
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 106
            movaps 3264(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3280(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 107
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 107
               movsldup 1584(pB), rB2    /* port 2 */
            #elif KB == 107
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 107
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 107
               movshdup 1584(pB), rB3    /* port 2 */
            #elif KB == 107
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 107
            movaps 3296(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3312(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 108
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 108
               movsldup 1600(pB), rB2    /* port 2 */
            #elif KB == 108
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 108
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 108
               movshdup 1600(pB), rB3    /* port 2 */
            #elif KB == 108
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 108
            movaps 3328(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3344(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 109
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 109
               movsldup 1616(pB), rB2    /* port 2 */
            #elif KB == 109
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 109
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 109
               movshdup 1616(pB), rB3    /* port 2 */
            #elif KB == 109
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 109
            movaps 3360(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3376(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 110
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 110
               movsldup 1632(pB), rB2    /* port 2 */
            #elif KB == 110
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 110
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 110
               movshdup 1632(pB), rB3    /* port 2 */
            #elif KB == 110
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 110
            movaps 3392(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3408(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 111
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 111
               movsldup 1648(pB), rB2    /* port 2 */
            #elif KB == 111
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 111
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 111
               movshdup 1648(pB), rB3    /* port 2 */
            #elif KB == 111
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 111
            movaps 3424(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3440(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 112
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 112
               movsldup 1664(pB), rB2    /* port 2 */
            #elif KB == 112
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 112
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 112
               movshdup 1664(pB), rB3    /* port 2 */
            #elif KB == 112
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 112
            movaps 3456(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3472(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 113
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 113
               movsldup 1680(pB), rB2    /* port 2 */
            #elif KB == 113
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 113
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 113
               movshdup 1680(pB), rB3    /* port 2 */
            #elif KB == 113
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 113
            movaps 3488(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3504(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 114
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 114
               movsldup 1696(pB), rB2    /* port 2 */
            #elif KB == 114
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 114
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 114
               movshdup 1696(pB), rB3    /* port 2 */
            #elif KB == 114
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 114
            movaps 3520(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3536(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 115
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 115
               movsldup 1712(pB), rB2    /* port 2 */
            #elif KB == 115
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 115
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 115
               movshdup 1712(pB), rB3    /* port 2 */
            #elif KB == 115
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 115
            movaps 3552(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3568(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 116
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 116
               movsldup 1728(pB), rB2    /* port 2 */
            #elif KB == 116
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 116
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 116
               movshdup 1728(pB), rB3    /* port 2 */
            #elif KB == 116
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 116
            movaps 3584(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3600(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 117
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 117
               movsldup 1744(pB), rB2    /* port 2 */
            #elif KB == 117
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 117
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 117
               movshdup 1744(pB), rB3    /* port 2 */
            #elif KB == 117
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 117
            movaps 3616(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3632(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 118
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 118
               movsldup 1760(pB), rB2    /* port 2 */
            #elif KB == 118
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 118
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 118
               movshdup 1760(pB), rB3    /* port 2 */
            #elif KB == 118
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 118
            movaps 3648(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3664(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 119
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 119
               movsldup 1776(pB), rB2    /* port 2 */
            #elif KB == 119
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 119
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 119
               movshdup 1776(pB), rB3    /* port 2 */
            #elif KB == 119
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 119
            movaps 3680(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3696(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 120
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 120
               movsldup 1792(pB), rB2    /* port 2 */
            #elif KB == 120
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 120
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 120
               movshdup 1792(pB), rB3    /* port 2 */
            #elif KB == 120
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 120
            movaps 3712(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3728(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 121
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 121
               movsldup 1808(pB), rB2    /* port 2 */
            #elif KB == 121
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 121
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 121
               movshdup 1808(pB), rB3    /* port 2 */
            #elif KB == 121
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 121
            movaps 3744(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3760(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 122
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 122
               movsldup 1824(pB), rB2    /* port 2 */
            #elif KB == 122
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 122
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 122
               movshdup 1824(pB), rB3    /* port 2 */
            #elif KB == 122
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 122
            movaps 3776(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3792(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 123
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 123
               movsldup 1840(pB), rB2    /* port 2 */
            #elif KB == 123
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 123
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 123
               movshdup 1840(pB), rB3    /* port 2 */
            #elif KB == 123
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 123
            movaps 3808(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3824(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 124
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 124
               movsldup 1856(pB), rB2    /* port 2 */
            #elif KB == 124
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 124
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 124
               movshdup 1856(pB), rB3    /* port 2 */
            #elif KB == 124
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 124
            movaps 3840(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3856(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 125
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 125
               movsldup 1872(pB), rB2    /* port 2 */
            #elif KB == 125
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 125
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 125
               movshdup 1872(pB), rB3    /* port 2 */
            #elif KB == 125
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 125
            movaps 3872(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3888(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 126
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 126
               movsldup 1888(pB), rB2    /* port 2 */
            #elif KB == 126
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 126
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 126
               movshdup 1888(pB), rB3    /* port 2 */
            #elif KB == 126
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 126
            movaps 3904(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3920(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 127
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 127
               movsldup 1904(pB), rB2    /* port 2 */
            #elif KB == 127
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 127
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 127
               movshdup 1904(pB), rB3    /* port 2 */
            #elif KB == 127
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 127
            movaps 3936(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3952(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 128
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 128
               movsldup 1920(pB), rB2    /* port 2 */
            #elif KB == 128
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 128
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 128
               movshdup 1920(pB), rB3    /* port 2 */
            #elif KB == 128
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 128
            movaps 3968(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3984(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 129
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 129
               movsldup 1936(pB), rB2    /* port 2 */
            #elif KB == 129
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 129
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 129
               movshdup 1936(pB), rB3    /* port 2 */
            #elif KB == 129
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 129
            movaps 4000(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4016(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 130
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 130
               movsldup 1952(pB), rB2    /* port 2 */
            #elif KB == 130
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 130
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 130
               movshdup 1952(pB), rB3    /* port 2 */
            #elif KB == 130
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 130
            movaps 4032(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4048(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 131
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 131
               movsldup 1968(pB), rB2    /* port 2 */
            #elif KB == 131
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 131
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 131
               movshdup 1968(pB), rB3    /* port 2 */
            #elif KB == 131
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 131
            movaps 4064(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4080(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 132
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 132
               movsldup 1984(pB), rB2    /* port 2 */
            #elif KB == 132
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 132
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 132
               movshdup 1984(pB), rB3    /* port 2 */
            #elif KB == 132
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 132
            movaps 4096(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4112(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 133
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 133
               movsldup 2000(pB), rB2    /* port 2 */
            #elif KB == 133
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 133
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 133
               movshdup 2000(pB), rB3    /* port 2 */
            #elif KB == 133
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 133
            movaps 4128(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4144(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 134
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 134
               movsldup 2016(pB), rB2    /* port 2 */
            #elif KB == 134
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 134
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 134
               movshdup 2016(pB), rB3    /* port 2 */
            #elif KB == 134
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 134
            movaps 4160(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4176(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 135
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 135
               movsldup 2032(pB), rB2    /* port 2 */
            #elif KB == 135
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 135
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 135
               movshdup 2032(pB), rB3    /* port 2 */
            #elif KB == 135
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 135
            movaps 4192(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4208(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 136
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 136
               movsldup 2048(pB), rB2    /* port 2 */
            #elif KB == 136
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 136
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 136
               movshdup 2048(pB), rB3    /* port 2 */
            #elif KB == 136
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 136
            movaps 4224(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4240(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 137
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 137
               movsldup 2064(pB), rB2    /* port 2 */
            #elif KB == 137
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 137
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 137
               movshdup 2064(pB), rB3    /* port 2 */
            #elif KB == 137
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 137
            movaps 4256(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4272(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 138
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 138
               movsldup 2080(pB), rB2    /* port 2 */
            #elif KB == 138
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 138
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 138
               movshdup 2080(pB), rB3    /* port 2 */
            #elif KB == 138
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 138
            movaps 4288(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4304(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 139
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 139
               movsldup 2096(pB), rB2    /* port 2 */
            #elif KB == 139
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 139
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 139
               movshdup 2096(pB), rB3    /* port 2 */
            #elif KB == 139
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 139
            movaps 4320(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4336(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 140
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 140
               movsldup 2112(pB), rB2    /* port 2 */
            #elif KB == 140
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 140
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 140
               movshdup 2112(pB), rB3    /* port 2 */
            #elif KB == 140
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 140
            movaps 4352(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4368(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 141
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 141
               movsldup 2128(pB), rB2    /* port 2 */
            #elif KB == 141
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 141
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 141
               movshdup 2128(pB), rB3    /* port 2 */
            #elif KB == 141
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 141
            movaps 4384(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4400(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 142
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 142
               movsldup 2144(pB), rB2    /* port 2 */
            #elif KB == 142
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 142
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 142
               movshdup 2144(pB), rB3    /* port 2 */
            #elif KB == 142
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 142
            movaps 4416(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4432(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 143
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 143
               movsldup 2160(pB), rB2    /* port 2 */
            #elif KB == 143
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 143
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 143
               movshdup 2160(pB), rB3    /* port 2 */
            #elif KB == 143
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 143
            movaps 4448(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4464(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 144
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 144
               movsldup 2176(pB), rB2    /* port 2 */
            #elif KB == 144
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 144
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 144
               movshdup 2176(pB), rB3    /* port 2 */
            #elif KB == 144
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 144
            movaps 4480(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4496(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 145
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 145
               movsldup 2192(pB), rB2    /* port 2 */
            #elif KB == 145
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 145
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 145
               movshdup 2192(pB), rB3    /* port 2 */
            #elif KB == 145
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 145
            movaps 4512(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4528(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 146
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 146
               movsldup 2208(pB), rB2    /* port 2 */
            #elif KB == 146
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 146
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 146
               movshdup 2208(pB), rB3    /* port 2 */
            #elif KB == 146
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 146
            movaps 4544(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4560(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 147
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 147
               movsldup 2224(pB), rB2    /* port 2 */
            #elif KB == 147
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 147
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 147
               movshdup 2224(pB), rB3    /* port 2 */
            #elif KB == 147
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 147
            movaps 4576(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4592(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 148
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 148
               movsldup 2240(pB), rB2    /* port 2 */
            #elif KB == 148
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 148
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 148
               movshdup 2240(pB), rB3    /* port 2 */
            #elif KB == 148
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 148
            movaps 4608(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4624(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 149
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 149
               movsldup 2256(pB), rB2    /* port 2 */
            #elif KB == 149
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 149
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 149
               movshdup 2256(pB), rB3    /* port 2 */
            #elif KB == 149
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 149
            movaps 4640(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4656(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 150
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 150
               movsldup 2272(pB), rB2    /* port 2 */
            #elif KB == 150
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 150
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 150
               movshdup 2272(pB), rB3    /* port 2 */
            #elif KB == 150
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 150
            movaps 4672(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4688(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 151
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 151
               movsldup 2288(pB), rB2    /* port 2 */
            #elif KB == 151
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 151
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 151
               movshdup 2288(pB), rB3    /* port 2 */
            #elif KB == 151
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 151
            movaps 4704(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4720(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 152
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 152
               movsldup 2304(pB), rB2    /* port 2 */
            #elif KB == 152
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 152
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 152
               movshdup 2304(pB), rB3    /* port 2 */
            #elif KB == 152
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 152
            movaps 4736(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4752(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 153
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 153
               movsldup 2320(pB), rB2    /* port 2 */
            #elif KB == 153
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 153
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 153
               movshdup 2320(pB), rB3    /* port 2 */
            #elif KB == 153
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 153
            movaps 4768(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4784(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 154
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 154
               movsldup 2336(pB), rB2    /* port 2 */
            #elif KB == 154
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 154
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 154
               movshdup 2336(pB), rB3    /* port 2 */
            #elif KB == 154
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 154
            movaps 4800(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4816(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 155
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 155
               movsldup 2352(pB), rB2    /* port 2 */
            #elif KB == 155
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 155
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 155
               movshdup 2352(pB), rB3    /* port 2 */
            #elif KB == 155
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 155
            movaps 4832(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4848(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 156
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 156
               movsldup 2368(pB), rB2    /* port 2 */
            #elif KB == 156
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 156
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 156
               movshdup 2368(pB), rB3    /* port 2 */
            #elif KB == 156
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 156
            movaps 4864(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4880(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 157
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 157
               movsldup 2384(pB), rB2    /* port 2 */
            #elif KB == 157
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 157
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 157
               movshdup 2384(pB), rB3    /* port 2 */
            #elif KB == 157
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 157
            movaps 4896(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4912(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 158
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 158
               movsldup 2400(pB), rB2    /* port 2 */
            #elif KB == 158
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 158
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 158
               movshdup 2400(pB), rB3    /* port 2 */
            #elif KB == 158
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 158
            movaps 4928(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4944(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 159
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 159
               movsldup 2416(pB), rB2    /* port 2 */
            #elif KB == 159
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 159
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 159
               movshdup 2416(pB), rB3    /* port 2 */
            #elif KB == 159
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 159
            movaps 4960(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4976(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 160
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 160
               movsldup 2432(pB), rB2    /* port 2 */
            #elif KB == 160
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 160
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 160
               movshdup 2432(pB), rB3    /* port 2 */
            #elif KB == 160
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 160
            movaps 4992(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5008(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 161
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 161
               movsldup 2448(pB), rB2    /* port 2 */
            #elif KB == 161
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 161
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 161
               movshdup 2448(pB), rB3    /* port 2 */
            #elif KB == 161
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 161
            movaps 5024(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5040(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 162
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 162
               movsldup 2464(pB), rB2    /* port 2 */
            #elif KB == 162
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 162
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 162
               movshdup 2464(pB), rB3    /* port 2 */
            #elif KB == 162
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 162
            movaps 5056(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5072(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 163
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 163
               movsldup 2480(pB), rB2    /* port 2 */
            #elif KB == 163
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 163
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 163
               movshdup 2480(pB), rB3    /* port 2 */
            #elif KB == 163
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 163
            movaps 5088(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5104(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 164
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 164
               movsldup 2496(pB), rB2    /* port 2 */
            #elif KB == 164
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 164
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 164
               movshdup 2496(pB), rB3    /* port 2 */
            #elif KB == 164
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 164
            movaps 5120(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5136(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 165
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 165
               movsldup 2512(pB), rB2    /* port 2 */
            #elif KB == 165
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 165
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 165
               movshdup 2512(pB), rB3    /* port 2 */
            #elif KB == 165
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 165
            movaps 5152(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5168(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 166
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 166
               movsldup 2528(pB), rB2    /* port 2 */
            #elif KB == 166
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 166
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 166
               movshdup 2528(pB), rB3    /* port 2 */
            #elif KB == 166
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 166
            movaps 5184(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5200(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 167
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 167
               movsldup 2544(pB), rB2    /* port 2 */
            #elif KB == 167
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 167
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 167
               movshdup 2544(pB), rB3    /* port 2 */
            #elif KB == 167
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 167
            movaps 5216(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5232(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 168
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 168
               movsldup 2560(pB), rB2    /* port 2 */
            #elif KB == 168
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 168
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 168
               movshdup 2560(pB), rB3    /* port 2 */
            #elif KB == 168
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 168
            movaps 5248(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5264(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 169
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 169
               movsldup 2576(pB), rB2    /* port 2 */
            #elif KB == 169
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 169
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 169
               movshdup 2576(pB), rB3    /* port 2 */
            #elif KB == 169
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 169
            movaps 5280(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5296(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 170
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 170
               movsldup 2592(pB), rB2    /* port 2 */
            #elif KB == 170
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 170
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 170
               movshdup 2592(pB), rB3    /* port 2 */
            #elif KB == 170
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 170
            movaps 5312(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5328(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 171
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 171
               movsldup 2608(pB), rB2    /* port 2 */
            #elif KB == 171
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 171
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 171
               movshdup 2608(pB), rB3    /* port 2 */
            #elif KB == 171
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 171
            movaps 5344(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5360(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 172
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 172
               movsldup 2624(pB), rB2    /* port 2 */
            #elif KB == 172
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 172
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 172
               movshdup 2624(pB), rB3    /* port 2 */
            #elif KB == 172
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 172
            movaps 5376(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5392(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 173
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 173
               movsldup 2640(pB), rB2    /* port 2 */
            #elif KB == 173
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 173
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 173
               movshdup 2640(pB), rB3    /* port 2 */
            #elif KB == 173
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 173
            movaps 5408(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5424(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 174
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 174
               movsldup 2656(pB), rB2    /* port 2 */
            #elif KB == 174
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 174
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 174
               movshdup 2656(pB), rB3    /* port 2 */
            #elif KB == 174
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 174
            movaps 5440(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5456(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 175
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 175
               movsldup 2672(pB), rB2    /* port 2 */
            #elif KB == 175
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 175
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 175
               movshdup 2672(pB), rB3    /* port 2 */
            #elif KB == 175
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 175
            movaps 5472(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5488(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 176
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 176
               movsldup 2688(pB), rB2    /* port 2 */
            #elif KB == 176
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 176
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 176
               movshdup 2688(pB), rB3    /* port 2 */
            #elif KB == 176
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 176
            movaps 5504(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5520(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 177
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 177
               movsldup 2704(pB), rB2    /* port 2 */
            #elif KB == 177
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 177
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 177
               movshdup 2704(pB), rB3    /* port 2 */
            #elif KB == 177
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 177
            movaps 5536(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5552(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 178
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 178
               movsldup 2720(pB), rB2    /* port 2 */
            #elif KB == 178
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 178
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 178
               movshdup 2720(pB), rB3    /* port 2 */
            #elif KB == 178
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 178
            movaps 5568(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5584(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 179
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 179
               movsldup 2736(pB), rB2    /* port 2 */
            #elif KB == 179
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 179
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 179
               movshdup 2736(pB), rB3    /* port 2 */
            #elif KB == 179
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 179
            movaps 5600(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5616(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 180
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 180
               movsldup 2752(pB), rB2    /* port 2 */
            #elif KB == 180
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 180
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 180
               movshdup 2752(pB), rB3    /* port 2 */
            #elif KB == 180
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 180
            movaps 5632(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5648(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 181
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 181
               movsldup 2768(pB), rB2    /* port 2 */
            #elif KB == 181
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 181
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 181
               movshdup 2768(pB), rB3    /* port 2 */
            #elif KB == 181
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 181
            movaps 5664(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5680(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 182
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 182
               movsldup 2784(pB), rB2    /* port 2 */
            #elif KB == 182
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 182
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 182
               movshdup 2784(pB), rB3    /* port 2 */
            #elif KB == 182
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 182
            movaps 5696(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5712(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 183
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 183
               movsldup 2800(pB), rB2    /* port 2 */
            #elif KB == 183
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 183
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 183
               movshdup 2800(pB), rB3    /* port 2 */
            #elif KB == 183
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 183
            movaps 5728(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5744(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 184
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 184
               movsldup 2816(pB), rB2    /* port 2 */
            #elif KB == 184
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 184
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 184
               movshdup 2816(pB), rB3    /* port 2 */
            #elif KB == 184
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 184
            movaps 5760(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5776(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 185
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 185
               movsldup 2832(pB), rB2    /* port 2 */
            #elif KB == 185
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 185
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 185
               movshdup 2832(pB), rB3    /* port 2 */
            #elif KB == 185
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 185
            movaps 5792(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5808(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 186
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 186
               movsldup 2848(pB), rB2    /* port 2 */
            #elif KB == 186
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 186
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 186
               movshdup 2848(pB), rB3    /* port 2 */
            #elif KB == 186
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 186
            movaps 5824(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5840(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 187
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 187
               movsldup 2864(pB), rB2    /* port 2 */
            #elif KB == 187
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 187
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 187
               movshdup 2864(pB), rB3    /* port 2 */
            #elif KB == 187
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 187
            movaps 5856(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5872(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 188
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 188
               movsldup 2880(pB), rB2    /* port 2 */
            #elif KB == 188
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 188
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 188
               movshdup 2880(pB), rB3    /* port 2 */
            #elif KB == 188
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 188
            movaps 5888(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5904(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 189
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 189
               movsldup 2896(pB), rB2    /* port 2 */
            #elif KB == 189
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 189
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 189
               movshdup 2896(pB), rB3    /* port 2 */
            #elif KB == 189
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 189
            movaps 5920(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5936(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 190
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 190
               movsldup 2912(pB), rB2    /* port 2 */
            #elif KB == 190
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 190
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 190
               movshdup 2912(pB), rB3    /* port 2 */
            #elif KB == 190
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 190
            movaps 5952(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5968(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 191
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 191
               movsldup 2928(pB), rB2    /* port 2 */
            #elif KB == 191
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 191
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 191
               movshdup 2928(pB), rB3    /* port 2 */
            #elif KB == 191
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 191
            movaps 5984(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6000(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 192
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 192
               movsldup 2944(pB), rB2    /* port 2 */
            #elif KB == 192
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 192
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 192
               movshdup 2944(pB), rB3    /* port 2 */
            #elif KB == 192
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 192
            movaps 6016(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6032(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 193
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 193
               movsldup 2960(pB), rB2    /* port 2 */
            #elif KB == 193
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 193
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 193
               movshdup 2960(pB), rB3    /* port 2 */
            #elif KB == 193
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 193
            movaps 6048(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6064(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 194
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 194
               movsldup 2976(pB), rB2    /* port 2 */
            #elif KB == 194
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 194
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 194
               movshdup 2976(pB), rB3    /* port 2 */
            #elif KB == 194
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 194
            movaps 6080(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6096(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 195
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 195
               movsldup 2992(pB), rB2    /* port 2 */
            #elif KB == 195
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 195
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 195
               movshdup 2992(pB), rB3    /* port 2 */
            #elif KB == 195
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 195
            movaps 6112(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6128(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 196
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 196
               movsldup 3008(pB), rB2    /* port 2 */
            #elif KB == 196
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 196
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 196
               movshdup 3008(pB), rB3    /* port 2 */
            #elif KB == 196
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 196
            movaps 6144(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6160(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 197
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 197
               movsldup 3024(pB), rB2    /* port 2 */
            #elif KB == 197
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 197
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 197
               movshdup 3024(pB), rB3    /* port 2 */
            #elif KB == 197
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 197
            movaps 6176(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6192(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 198
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 198
               movsldup 3040(pB), rB2    /* port 2 */
            #elif KB == 198
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 198
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 198
               movshdup 3040(pB), rB3    /* port 2 */
            #elif KB == 198
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 198
            movaps 6208(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6224(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 199
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 199
               movsldup 3056(pB), rB2    /* port 2 */
            #elif KB == 199
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 199
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 199
               movshdup 3056(pB), rB3    /* port 2 */
            #elif KB == 199
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 199
            movaps 6240(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6256(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 200
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 200
               movsldup 3072(pB), rB2    /* port 2 */
            #elif KB == 200
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 200
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 200
               movshdup 3072(pB), rB3    /* port 2 */
            #elif KB == 200
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 200
            movaps 6272(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6288(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 201
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 201
               movsldup 3088(pB), rB2    /* port 2 */
            #elif KB == 201
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 201
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 201
               movshdup 3088(pB), rB3    /* port 2 */
            #elif KB == 201
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 201
            movaps 6304(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6320(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 202
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 202
               movsldup 3104(pB), rB2    /* port 2 */
            #elif KB == 202
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 202
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 202
               movshdup 3104(pB), rB3    /* port 2 */
            #elif KB == 202
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 202
            movaps 6336(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6352(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 203
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 203
               movsldup 3120(pB), rB2    /* port 2 */
            #elif KB == 203
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 203
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 203
               movshdup 3120(pB), rB3    /* port 2 */
            #elif KB == 203
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 203
            movaps 6368(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6384(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 204
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 204
               movsldup 3136(pB), rB2    /* port 2 */
            #elif KB == 204
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 204
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 204
               movshdup 3136(pB), rB3    /* port 2 */
            #elif KB == 204
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 204
            movaps 6400(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6416(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 205
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 205
               movsldup 3152(pB), rB2    /* port 2 */
            #elif KB == 205
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 205
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 205
               movshdup 3152(pB), rB3    /* port 2 */
            #elif KB == 205
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 205
            movaps 6432(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6448(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 206
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 206
               movsldup 3168(pB), rB2    /* port 2 */
            #elif KB == 206
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 206
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 206
               movshdup 3168(pB), rB3    /* port 2 */
            #elif KB == 206
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 206
            movaps 6464(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6480(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 207
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 207
               movsldup 3184(pB), rB2    /* port 2 */
            #elif KB == 207
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 207
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 207
               movshdup 3184(pB), rB3    /* port 2 */
            #elif KB == 207
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 207
            movaps 6496(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6512(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 208
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 208
               movsldup 3200(pB), rB2    /* port 2 */
            #elif KB == 208
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 208
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 208
               movshdup 3200(pB), rB3    /* port 2 */
            #elif KB == 208
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 208
            movaps 6528(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6544(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 209
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 209
               movsldup 3216(pB), rB2    /* port 2 */
            #elif KB == 209
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 209
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 209
               movshdup 3216(pB), rB3    /* port 2 */
            #elif KB == 209
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 209
            movaps 6560(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6576(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 210
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 210
               movsldup 3232(pB), rB2    /* port 2 */
            #elif KB == 210
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 210
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 210
               movshdup 3232(pB), rB3    /* port 2 */
            #elif KB == 210
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 210
            movaps 6592(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6608(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 211
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 211
               movsldup 3248(pB), rB2    /* port 2 */
            #elif KB == 211
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 211
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 211
               movshdup 3248(pB), rB3    /* port 2 */
            #elif KB == 211
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 211
            movaps 6624(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6640(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 212
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 212
               movsldup 3264(pB), rB2    /* port 2 */
            #elif KB == 212
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 212
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 212
               movshdup 3264(pB), rB3    /* port 2 */
            #elif KB == 212
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 212
            movaps 6656(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6672(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 213
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 213
               movsldup 3280(pB), rB2    /* port 2 */
            #elif KB == 213
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 213
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 213
               movshdup 3280(pB), rB3    /* port 2 */
            #elif KB == 213
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 213
            movaps 6688(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6704(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 214
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 214
               movsldup 3296(pB), rB2    /* port 2 */
            #elif KB == 214
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 214
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 214
               movshdup 3296(pB), rB3    /* port 2 */
            #elif KB == 214
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 214
            movaps 6720(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6736(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 215
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 215
               movsldup 3312(pB), rB2    /* port 2 */
            #elif KB == 215
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 215
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 215
               movshdup 3312(pB), rB3    /* port 2 */
            #elif KB == 215
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 215
            movaps 6752(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6768(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 216
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 216
               movsldup 3328(pB), rB2    /* port 2 */
            #elif KB == 216
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 216
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 216
               movshdup 3328(pB), rB3    /* port 2 */
            #elif KB == 216
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 216
            movaps 6784(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6800(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 217
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 217
               movsldup 3344(pB), rB2    /* port 2 */
            #elif KB == 217
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 217
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 217
               movshdup 3344(pB), rB3    /* port 2 */
            #elif KB == 217
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 217
            movaps 6816(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6832(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 218
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 218
               movsldup 3360(pB), rB2    /* port 2 */
            #elif KB == 218
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 218
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 218
               movshdup 3360(pB), rB3    /* port 2 */
            #elif KB == 218
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 218
            movaps 6848(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6864(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 219
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 219
               movsldup 3376(pB), rB2    /* port 2 */
            #elif KB == 219
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 219
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 219
               movshdup 3376(pB), rB3    /* port 2 */
            #elif KB == 219
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 219
            movaps 6880(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6896(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 220
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 220
               movsldup 3392(pB), rB2    /* port 2 */
            #elif KB == 220
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 220
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 220
               movshdup 3392(pB), rB3    /* port 2 */
            #elif KB == 220
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 220
            movaps 6912(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6928(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 221
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 221
               movsldup 3408(pB), rB2    /* port 2 */
            #elif KB == 221
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 221
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 221
               movshdup 3408(pB), rB3    /* port 2 */
            #elif KB == 221
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 221
            movaps 6944(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6960(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 222
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 222
               movsldup 3424(pB), rB2    /* port 2 */
            #elif KB == 222
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 222
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 222
               movshdup 3424(pB), rB3    /* port 2 */
            #elif KB == 222
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 222
            movaps 6976(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6992(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 223
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 223
               movsldup 3440(pB), rB2    /* port 2 */
            #elif KB == 223
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 223
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 223
               movshdup 3440(pB), rB3    /* port 2 */
            #elif KB == 223
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 223
            movaps 7008(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7024(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 224
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 224
               movsldup 3456(pB), rB2    /* port 2 */
            #elif KB == 224
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 224
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 224
               movshdup 3456(pB), rB3    /* port 2 */
            #elif KB == 224
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 224
            movaps 7040(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7056(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 225
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 225
               movsldup 3472(pB), rB2    /* port 2 */
            #elif KB == 225
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 225
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 225
               movshdup 3472(pB), rB3    /* port 2 */
            #elif KB == 225
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 225
            movaps 7072(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7088(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 226
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 226
               movsldup 3488(pB), rB2    /* port 2 */
            #elif KB == 226
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 226
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 226
               movshdup 3488(pB), rB3    /* port 2 */
            #elif KB == 226
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 226
            movaps 7104(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7120(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 227
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 227
               movsldup 3504(pB), rB2    /* port 2 */
            #elif KB == 227
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 227
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 227
               movshdup 3504(pB), rB3    /* port 2 */
            #elif KB == 227
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 227
            movaps 7136(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7152(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 228
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 228
               movsldup 3520(pB), rB2    /* port 2 */
            #elif KB == 228
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 228
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 228
               movshdup 3520(pB), rB3    /* port 2 */
            #elif KB == 228
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 228
            movaps 7168(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7184(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 229
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 229
               movsldup 3536(pB), rB2    /* port 2 */
            #elif KB == 229
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 229
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 229
               movshdup 3536(pB), rB3    /* port 2 */
            #elif KB == 229
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 229
            movaps 7200(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7216(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 230
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 230
               movsldup 3552(pB), rB2    /* port 2 */
            #elif KB == 230
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 230
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 230
               movshdup 3552(pB), rB3    /* port 2 */
            #elif KB == 230
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 230
            movaps 7232(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7248(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 231
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 231
               movsldup 3568(pB), rB2    /* port 2 */
            #elif KB == 231
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 231
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 231
               movshdup 3568(pB), rB3    /* port 2 */
            #elif KB == 231
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 231
            movaps 7264(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7280(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 232
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 232
               movsldup 3584(pB), rB2    /* port 2 */
            #elif KB == 232
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 232
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 232
               movshdup 3584(pB), rB3    /* port 2 */
            #elif KB == 232
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 232
            movaps 7296(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7312(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 233
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 233
               movsldup 3600(pB), rB2    /* port 2 */
            #elif KB == 233
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 233
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 233
               movshdup 3600(pB), rB3    /* port 2 */
            #elif KB == 233
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 233
            movaps 7328(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7344(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 234
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 234
               movsldup 3616(pB), rB2    /* port 2 */
            #elif KB == 234
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 234
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 234
               movshdup 3616(pB), rB3    /* port 2 */
            #elif KB == 234
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 234
            movaps 7360(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7376(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 235
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 235
               movsldup 3632(pB), rB2    /* port 2 */
            #elif KB == 235
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 235
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 235
               movshdup 3632(pB), rB3    /* port 2 */
            #elif KB == 235
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 235
            movaps 7392(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7408(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 236
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 236
               movsldup 3648(pB), rB2    /* port 2 */
            #elif KB == 236
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 236
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 236
               movshdup 3648(pB), rB3    /* port 2 */
            #elif KB == 236
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 236
            movaps 7424(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7440(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 237
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 237
               movsldup 3664(pB), rB2    /* port 2 */
            #elif KB == 237
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 237
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 237
               movshdup 3664(pB), rB3    /* port 2 */
            #elif KB == 237
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 237
            movaps 7456(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7472(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 238
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 238
               movsldup 3680(pB), rB2    /* port 2 */
            #elif KB == 238
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 238
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 238
               movshdup 3680(pB), rB3    /* port 2 */
            #elif KB == 238
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 238
            movaps 7488(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7504(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 239
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 239
               movsldup 3696(pB), rB2    /* port 2 */
            #elif KB == 239
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 239
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 239
               movshdup 3696(pB), rB3    /* port 2 */
            #elif KB == 239
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 239
            movaps 7520(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7536(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 240
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 240
               movsldup 3712(pB), rB2    /* port 2 */
            #elif KB == 240
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 240
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 240
               movshdup 3712(pB), rB3    /* port 2 */
            #elif KB == 240
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 240
            movaps 7552(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7568(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 241
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 241
               movsldup 3728(pB), rB2    /* port 2 */
            #elif KB == 241
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 241
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 241
               movshdup 3728(pB), rB3    /* port 2 */
            #elif KB == 241
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 241
            movaps 7584(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7600(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 242
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 242
               movsldup 3744(pB), rB2    /* port 2 */
            #elif KB == 242
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 242
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 242
               movshdup 3744(pB), rB3    /* port 2 */
            #elif KB == 242
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 242
            movaps 7616(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7632(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 243
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 243
               movsldup 3760(pB), rB2    /* port 2 */
            #elif KB == 243
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 243
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 243
               movshdup 3760(pB), rB3    /* port 2 */
            #elif KB == 243
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 243
            movaps 7648(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7664(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 244
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 244
               movsldup 3776(pB), rB2    /* port 2 */
            #elif KB == 244
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 244
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 244
               movshdup 3776(pB), rB3    /* port 2 */
            #elif KB == 244
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 244
            movaps 7680(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7696(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 245
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 245
               movsldup 3792(pB), rB2    /* port 2 */
            #elif KB == 245
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 245
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 245
               movshdup 3792(pB), rB3    /* port 2 */
            #elif KB == 245
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 245
            movaps 7712(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7728(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 246
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 246
               movsldup 3808(pB), rB2    /* port 2 */
            #elif KB == 246
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 246
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 246
               movshdup 3808(pB), rB3    /* port 2 */
            #elif KB == 246
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 246
            movaps 7744(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7760(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 247
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 247
               movsldup 3824(pB), rB2    /* port 2 */
            #elif KB == 247
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 247
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 247
               movshdup 3824(pB), rB3    /* port 2 */
            #elif KB == 247
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 247
            movaps 7776(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7792(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 248
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 248
               movsldup 3840(pB), rB2    /* port 2 */
            #elif KB == 248
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 248
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 248
               movshdup 3840(pB), rB3    /* port 2 */
            #elif KB == 248
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 248
            movaps 7808(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7824(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 249
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 249
               movsldup 3856(pB), rB2    /* port 2 */
            #elif KB == 249
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 249
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 249
               movshdup 3856(pB), rB3    /* port 2 */
            #elif KB == 249
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 249
            movaps 7840(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7856(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 250
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 250
               movsldup 3872(pB), rB2    /* port 2 */
            #elif KB == 250
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 250
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 250
               movshdup 3872(pB), rB3    /* port 2 */
            #elif KB == 250
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 250
            movaps 7872(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7888(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 251
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 251
               movsldup 3888(pB), rB2    /* port 2 */
            #elif KB == 251
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 251
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 251
               movshdup 3888(pB), rB3    /* port 2 */
            #elif KB == 251
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 251
            movaps 7904(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7920(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 252
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 252
               movsldup 3904(pB), rB2    /* port 2 */
            #elif KB == 252
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 252
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 252
               movshdup 3904(pB), rB3    /* port 2 */
            #elif KB == 252
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 252
            movaps 7936(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7952(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 253
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 253
               movsldup 3920(pB), rB2    /* port 2 */
            #elif KB == 253
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 253
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 253
               movshdup 3920(pB), rB3    /* port 2 */
            #elif KB == 253
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 253
            movaps 7968(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7984(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 254
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 254
               movsldup 3936(pB), rB2    /* port 2 */
            #elif KB == 254
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 254
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 254
               movshdup 3936(pB), rB3    /* port 2 */
            #elif KB == 254
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 254
            movaps 8000(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8016(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 255
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 255
               movsldup 3952(pB), rB2    /* port 2 */
            #elif KB == 255
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 255
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 255
               movshdup 3952(pB), rB3    /* port 2 */
            #elif KB == 255
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if KB > 255
            movaps 8032(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8048(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            movhlps rB2, rB2            /* port 5, {b2, b2, b2, b2} */
            mulps  rB2, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC02            /* port 1, 4 bytes 12 */
            #if KB == 256
               pref (pf)                /* port 2 */
            #endif

            mulps  rA1, rB2             /* port 0 */
            addps  rB2, rC12            /* port 1 */
            #if KB > 256
               movsldup 3968(pB), rB2    /* port 2 */
            #elif KB == 256
               prefB 64(pfB)            /* port 2 */
            #endif
            movaps rA0, rm0             /* port 5 */

            movddup rB3, rB1            /* port 5 */
            mulps  rB1, rm0             /* port 0 */
            addps  rm0, rC01            /* port 1 */

            movhlps rB3, rB3            /* port 5 */
            mulps  rA1, rB1             /* port 0 */
            addps  rB1, rC11            /* port 1 */

            #if KB > 256
               movddup rB2, rB0         /* port 5, {b0, b0, b0, b0} */
            #endif
            mulps rB3, rA0              /* port 0 */
            addps rA0, rC03             /* port 1 */

            mulps rB3, rA1              /* port 0 */
            addps rA1, rC13             /* port 1 */
            #if KB > 256
               movshdup 3968(pB), rB3    /* port 2 */
            #elif KB == 256
               add incPF, pf            /* port5 */
               add incPF, pfB           /* port0/1/5 */
            #endif
         #endif
         #if defined(BETA1) || defined(BETAN1)
            BETCOP (pC), rC00
            movaps rC00, (pC)
            BETCOP 16(pC), rC10
            movaps rC10, 16(pC)
            BETCOP 32(pC), rC01
            movaps rC01, 32(pC)
            BETCOP 48(pC), rC11
            movaps rC11, 48(pC)
            BETCOP 64(pC), rC02
            movaps rC02, 64(pC)
            BETCOP 80(pC), rC12
            movaps rC12, 80(pC)
            BETCOP 96(pC), rC03
            movaps rC03, 96(pC)
            BETCOP 112(pC), rC13
            movaps rC13, 112(pC)
         #else
            movaps rC00, (pC)
            movaps rC10, 16(pC)
            movaps rC01, 32(pC)
            movaps rC11, 48(pC)
            movaps rC02, 64(pC)
            movaps rC12, 80(pC)
            movaps rC03, 96(pC)
            movaps rC13, 112(pC)
         #endif
         sub $-128, pC
         add $KB*4*4, pB
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
   add $FSIZE, %rsp
   ret
