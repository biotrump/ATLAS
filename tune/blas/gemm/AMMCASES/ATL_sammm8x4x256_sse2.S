/*
 *             Automatically Tuned Linear Algebra Software v3.11.31
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
#define movaps movaps
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
#define pfC     %r13
#define r256    %r14

#define rm0     %xmm0
#define rA0     %xmm1
#define rA1     %xmm2
#define rB0     %xmm3
#define rB1     %xmm4
#define rB2     %xmm5
#define rC00    %xmm6
#define rC10    %xmm7
#define rC01    %xmm8
#define rC11    %xmm9
#define rC02    %xmm10
#define rC12    %xmm11
#define rC03    %xmm12
#define rC13    %xmm13
#ifndef pref
   #define pref prefetcht1
#endif
#ifndef prefB
   #define prefB prefetcht1
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
   movq FSIZE+16(%rsp), pfB     /* pf = pBn */
   movq FSIZE+8(%rsp), pfA      /* pfB = pAn */
   movq FSIZE+24(%rsp), pfC    /* pfC = pCn */
   mov $4*4*4, incPF
/*
 * Make it so pA/pB always use 32-bit constants by subtracting 128
 */
   sub $96, pA
   sub $96, pB
   mov $KB*8*4, incAm           /* incAm = KB*MU*size */
   movq pB, pB0

           ALIGN32
   .local MNLOOP
   MNLOOP:
/*
 *       Peel first iteration of K-loop to handle init of C to 0
 */
#if 0
         xorps rC00, rC00
         movaps rC00, rC10
         xorps rC01, rC01
         movaps rC00, rC11
         xorps rC02, rC02
         movaps rC00, rC12
         xorps rC03, rC03
         movaps rC00, rC13
#endif
         .byte 0x3e                     /* no prt, 1 bytes  1 */
         movaps 96(pB), rB2             /* port 2, 4 bytes  5 */
         .byte 0x3e                     /* no prt, 1 bytes  6 */
         pshufd $0x00, rB2, rB0         /* port 5, 5 bytes 11 */

         .byte 0x3e                     /* no prt, 1 bytes 12 */
         movaps 96(pA), rC00            /* port 2, 4 bytes 16 */
         movaps rC00, rC01              /* port 5, 4 bytes 04 */
         mulps  rB0, rC00               /* port 0, 3 bytes 07 */

         movaps 112(pA), rC10           /* port 2, 4 bytes 11 */
         .byte 0x3e                     /* no prt, 1 bytes 12 */
         movaps rC10, rC11              /* port 5, 4 bytes 16 */
         mulps rB0, rC10                /* port 0, 3 bytes 03 */

         pshufd $0x55, rB2, rB1         /* port 5, 5 bytes 08 */
            prefetcht0 (pC)             /* port 2, 4 bytes 12 */

         movaps rC01, rC02              /* port 5, 4 bytes 16 */
         mulps rB1, rC01                /* port 0, 4 bytes 04 */

         movaps rC11, rC12              /* port 5, 4 bytes 08 */
         mulps rB1, rC11                /* port 0, 4 bytes 12 */
            movaps 112(pB), rB1         /* port 2, 4 bytes 16 */

        movaps rC02, rC03               /* port 5, 4 bytes 04 */
         .byte 0x3e                     /* no prt, 1 bytes 05 */
           prefetcht0 64(pC)            /* port 2, 5 bytes 10 */

         .byte 0x3e                     /* no prt, 1 bytes 11 */
        pshufd $0xAA, rB2, rB0          /* port 5, 5 bytes 16 */
        mulps rB0, rC02                 /* port 0, 4 bytes 04 */
           prefetcht2 (pfB)             /* port 2, 3 bytes 07 */

        movaps rC12, rC13               /* port 5, 4 bytes 11 */
         .byte 0x3e                     /* no prt, 1 bytes 12 */
        mulps rB0, rC12                 /* port 0, 4 bytes 16 */
         .byte 0x3e                     /* no prt, 1 bytes 01 */
           prefetcht2 (pfA)             /* port 2, 4 bytes 05 */

         .byte 0x3e                     /* no prt, 1 bytes 06 */
        pshufd $0xFF, rB2, rm0          /* port 5, 5 bytes 11 */
         .byte 0x3e                     /* no prt, 1 bytes 12 */
        mulps rm0, rC03                 /* port 0, 4 bytes 16 */

            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */
        mulps rm0, rC13                 /* port 0, 4 bytes 09 */

            movaps rB1, rB2             /* port 5, 3 bytes 12 */
            .byte 0x66                  /* no prt, 1 bytes 13 */
            .byte 0x66                  /* no prt, 1 bytes 14 */
            .byte 0x66                  /* no prt, 1 bytes 15 */
            nop                         /* port 0, 1 bytes 16 */
/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 1
            movaps 128(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 144(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 128(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 2
            movaps 160(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 176(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 144(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 3
            movaps 192(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 208(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 160(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 4
            movaps 224(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 240(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 176(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 5
            movaps 256(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 272(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 192(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 6
            movaps 288(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 304(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 208(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 7
            movaps 320(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 336(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 224(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 8
            movaps 352(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 368(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 240(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 9
            movaps 384(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 400(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 256(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 10
            movaps 416(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 432(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 272(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 11
            movaps 448(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 464(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 288(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 12
            movaps 480(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 496(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 304(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 13
            movaps 512(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 528(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 320(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 14
            movaps 544(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 560(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 336(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 15
            movaps 576(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 592(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 352(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 16
            movaps 608(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 624(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 368(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 17
            movaps 640(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 656(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 384(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 18
            movaps 672(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 688(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 400(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 19
            movaps 704(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 720(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 416(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 20
            movaps 736(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 752(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 432(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 21
            movaps 768(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 784(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 448(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 22
            movaps 800(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 816(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 464(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 23
            movaps 832(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 848(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 480(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 24
            movaps 864(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 880(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 496(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 25
            movaps 896(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 912(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 512(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 26
            movaps 928(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 944(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 528(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 27
            movaps 960(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 976(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 544(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 28
            movaps 992(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1008(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 560(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 29
            movaps 1024(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1040(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 576(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 30
            movaps 1056(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1072(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 592(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 31
            movaps 1088(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1104(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 608(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 32
            movaps 1120(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1136(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 624(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 33
            movaps 1152(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1168(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 640(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 34
            movaps 1184(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1200(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 656(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 35
            movaps 1216(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1232(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 672(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 36
            movaps 1248(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1264(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 688(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 37
            movaps 1280(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1296(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 704(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 38
            movaps 1312(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1328(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 720(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 39
            movaps 1344(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1360(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 736(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 40
            movaps 1376(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1392(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 752(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 41
            movaps 1408(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1424(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 768(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 42
            movaps 1440(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1456(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 784(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 43
            movaps 1472(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1488(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 800(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 44
            movaps 1504(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1520(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 816(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 45
            movaps 1536(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1552(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 832(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 46
            movaps 1568(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1584(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 848(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 47
            movaps 1600(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1616(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 864(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 48
            movaps 1632(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1648(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 880(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 49
            movaps 1664(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1680(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 896(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 50
            movaps 1696(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1712(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 912(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 51
            movaps 1728(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1744(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 928(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 52
            movaps 1760(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1776(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 944(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 53
            movaps 1792(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1808(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 960(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 54
            movaps 1824(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1840(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 976(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 55
            movaps 1856(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1872(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 992(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 56
            movaps 1888(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1904(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1008(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 57
            movaps 1920(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1936(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1024(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 58
            movaps 1952(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 1968(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1040(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 59
            movaps 1984(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2000(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1056(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 60
            movaps 2016(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2032(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1072(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 61
            movaps 2048(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2064(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1088(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 62
            movaps 2080(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2096(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1104(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 63
            movaps 2112(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2128(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1120(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 64
            movaps 2144(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2160(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1136(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 65
            movaps 2176(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2192(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1152(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 66
            movaps 2208(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2224(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1168(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 67
            movaps 2240(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2256(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1184(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 68
            movaps 2272(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2288(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1200(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 69
            movaps 2304(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2320(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1216(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 70
            movaps 2336(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2352(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1232(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 71
            movaps 2368(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2384(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1248(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 72
            movaps 2400(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2416(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1264(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 73
            movaps 2432(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2448(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1280(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 74
            movaps 2464(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2480(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1296(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 75
            movaps 2496(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2512(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1312(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 76
            movaps 2528(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2544(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1328(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 77
            movaps 2560(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2576(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1344(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 78
            movaps 2592(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2608(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1360(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 79
            movaps 2624(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2640(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1376(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 80
            movaps 2656(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2672(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1392(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 81
            movaps 2688(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2704(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1408(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 82
            movaps 2720(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2736(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1424(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 83
            movaps 2752(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2768(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1440(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 84
            movaps 2784(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2800(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1456(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 85
            movaps 2816(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2832(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1472(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 86
            movaps 2848(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2864(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1488(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 87
            movaps 2880(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2896(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1504(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 88
            movaps 2912(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2928(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1520(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 89
            movaps 2944(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2960(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1536(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 90
            movaps 2976(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 2992(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1552(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 91
            movaps 3008(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3024(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1568(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 92
            movaps 3040(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3056(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1584(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 93
            movaps 3072(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3088(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1600(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 94
            movaps 3104(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3120(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1616(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 95
            movaps 3136(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3152(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1632(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 96
            movaps 3168(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3184(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1648(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 97
            movaps 3200(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3216(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1664(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 98
            movaps 3232(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3248(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1680(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 99
            movaps 3264(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3280(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1696(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 100
            movaps 3296(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3312(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1712(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 101
            movaps 3328(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3344(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1728(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 102
            movaps 3360(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3376(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1744(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 103
            movaps 3392(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3408(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1760(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 104
            movaps 3424(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3440(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1776(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 105
            movaps 3456(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3472(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1792(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 106
            movaps 3488(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3504(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1808(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 107
            movaps 3520(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3536(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1824(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 108
            movaps 3552(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3568(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1840(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 109
            movaps 3584(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3600(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1856(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 110
            movaps 3616(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3632(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1872(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 111
            movaps 3648(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3664(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1888(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 112
            movaps 3680(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3696(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1904(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 113
            movaps 3712(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3728(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1920(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 114
            movaps 3744(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3760(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1936(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 115
            movaps 3776(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3792(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1952(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 116
            movaps 3808(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3824(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1968(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 117
            movaps 3840(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3856(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 1984(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 118
            movaps 3872(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3888(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2000(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 119
            movaps 3904(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3920(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2016(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 120
            movaps 3936(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3952(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2032(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 121
            movaps 3968(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 3984(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2048(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 122
            movaps 4000(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4016(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2064(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 123
            movaps 4032(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4048(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2080(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 124
            movaps 4064(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4080(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2096(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 125
            movaps 4096(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4112(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2112(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 126
            movaps 4128(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4144(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2128(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 127
            movaps 4160(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4176(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2144(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 128
            movaps 4192(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4208(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2160(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 129
            movaps 4224(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4240(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2176(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 130
            movaps 4256(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4272(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2192(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 131
            movaps 4288(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4304(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2208(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 132
            movaps 4320(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4336(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2224(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 133
            movaps 4352(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4368(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2240(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 134
            movaps 4384(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4400(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2256(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 135
            movaps 4416(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4432(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2272(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 136
            movaps 4448(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4464(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2288(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 137
            movaps 4480(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4496(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2304(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 138
            movaps 4512(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4528(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2320(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 139
            movaps 4544(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4560(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2336(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 140
            movaps 4576(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4592(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2352(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 141
            movaps 4608(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4624(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2368(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 142
            movaps 4640(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4656(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2384(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 143
            movaps 4672(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4688(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2400(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 144
            movaps 4704(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4720(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2416(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 145
            movaps 4736(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4752(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2432(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 146
            movaps 4768(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4784(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2448(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 147
            movaps 4800(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4816(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2464(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 148
            movaps 4832(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4848(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2480(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 149
            movaps 4864(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4880(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2496(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 150
            movaps 4896(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4912(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2512(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 151
            movaps 4928(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4944(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2528(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 152
            movaps 4960(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 4976(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2544(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 153
            movaps 4992(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5008(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2560(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 154
            movaps 5024(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5040(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2576(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 155
            movaps 5056(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5072(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2592(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 156
            movaps 5088(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5104(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2608(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 157
            movaps 5120(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5136(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2624(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 158
            movaps 5152(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5168(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2640(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 159
            movaps 5184(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5200(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2656(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 160
            movaps 5216(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5232(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2672(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 161
            movaps 5248(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5264(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2688(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 162
            movaps 5280(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5296(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2704(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 163
            movaps 5312(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5328(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2720(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 164
            movaps 5344(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5360(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2736(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 165
            movaps 5376(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5392(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2752(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 166
            movaps 5408(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5424(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2768(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 167
            movaps 5440(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5456(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2784(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 168
            movaps 5472(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5488(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2800(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 169
            movaps 5504(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5520(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2816(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 170
            movaps 5536(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5552(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2832(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 171
            movaps 5568(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5584(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2848(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 172
            movaps 5600(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5616(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2864(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 173
            movaps 5632(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5648(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2880(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 174
            movaps 5664(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5680(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2896(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 175
            movaps 5696(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5712(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2912(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 176
            movaps 5728(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5744(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2928(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 177
            movaps 5760(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5776(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2944(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 178
            movaps 5792(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5808(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2960(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 179
            movaps 5824(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5840(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2976(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 180
            movaps 5856(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5872(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 2992(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 181
            movaps 5888(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5904(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3008(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 182
            movaps 5920(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5936(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3024(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 183
            movaps 5952(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 5968(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3040(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 184
            movaps 5984(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6000(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3056(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 185
            movaps 6016(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6032(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3072(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 186
            movaps 6048(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6064(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3088(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 187
            movaps 6080(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6096(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3104(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 188
            movaps 6112(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6128(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3120(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 189
            movaps 6144(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6160(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3136(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 190
            movaps 6176(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6192(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3152(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 191
            movaps 6208(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6224(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3168(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 192
            movaps 6240(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6256(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3184(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 193
            movaps 6272(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6288(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3200(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 194
            movaps 6304(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6320(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3216(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 195
            movaps 6336(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6352(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3232(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 196
            movaps 6368(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6384(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3248(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 197
            movaps 6400(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6416(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3264(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 198
            movaps 6432(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6448(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3280(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 199
            movaps 6464(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6480(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3296(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 200
            movaps 6496(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6512(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3312(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 201
            movaps 6528(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6544(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3328(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 202
            movaps 6560(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6576(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3344(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 203
            movaps 6592(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6608(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3360(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 204
            movaps 6624(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6640(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3376(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 205
            movaps 6656(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6672(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3392(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 206
            movaps 6688(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6704(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3408(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 207
            movaps 6720(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6736(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3424(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 208
            movaps 6752(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6768(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3440(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 209
            movaps 6784(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6800(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3456(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 210
            movaps 6816(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6832(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3472(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 211
            movaps 6848(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6864(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3488(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 212
            movaps 6880(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6896(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3504(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 213
            movaps 6912(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6928(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3520(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 214
            movaps 6944(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6960(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3536(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 215
            movaps 6976(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 6992(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3552(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 216
            movaps 7008(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7024(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3568(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 217
            movaps 7040(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7056(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3584(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 218
            movaps 7072(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7088(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3600(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 219
            movaps 7104(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7120(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3616(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 220
            movaps 7136(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7152(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3632(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 221
            movaps 7168(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7184(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3648(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 222
            movaps 7200(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7216(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3664(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 223
            movaps 7232(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7248(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3680(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 224
            movaps 7264(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7280(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3696(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 225
            movaps 7296(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7312(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3712(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 226
            movaps 7328(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7344(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3728(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 227
            movaps 7360(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7376(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3744(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 228
            movaps 7392(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7408(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3760(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 229
            movaps 7424(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7440(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3776(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 230
            movaps 7456(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7472(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3792(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 231
            movaps 7488(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7504(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3808(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 232
            movaps 7520(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7536(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3824(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 233
            movaps 7552(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7568(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3840(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 234
            movaps 7584(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7600(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3856(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 235
            movaps 7616(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7632(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3872(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 236
            movaps 7648(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7664(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3888(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 237
            movaps 7680(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7696(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3904(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 238
            movaps 7712(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7728(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3920(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 239
            movaps 7744(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7760(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3936(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 240
            movaps 7776(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7792(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3952(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 241
            movaps 7808(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7824(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3968(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 242
            movaps 7840(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7856(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 3984(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 243
            movaps 7872(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7888(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4000(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 244
            movaps 7904(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7920(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4016(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 245
            movaps 7936(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7952(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4032(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 246
            movaps 7968(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 7984(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4048(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 247
            movaps 8000(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8016(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4064(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 248
            movaps 8032(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8048(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4080(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 249
            movaps 8064(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8080(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4096(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 250
            movaps 8096(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8112(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4112(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 251
            movaps 8128(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8144(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4128(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 252
            movaps 8160(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8176(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4144(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 253
            movaps 8192(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8208(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4160(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 254
            movaps 8224(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8240(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4176(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
         #endif
         #if KB > 255
            movaps 8256(pA), rA0         /* port 2, 7 bytes  7 */
            movaps rA0, rm0             /* port 5, 3 bytes 10 */
            mulps  rB0, rm0             /* port 0, 3 bytes 13 */
            addps  rm0, rC00            /* port 1, 3 bytes 16 */

            movaps 8272(pA), rA1         /* port 2, 7 bytes  7 */
            mulps  rA1, rB0             /* port 0, 3 bytes 10 */
            addps  rB0, rC10            /* port 1, 3 bytes 13 */
            movaps rA0, rm0             /* port 5, 3 bytes 16 */

            pshufd $0x55, rB2, rB1      /* port 5, 5 bytes  5 */
            mulps  rB1, rm0             /* port 0, 3 bytes  8 */
            addps  rm0, rC01            /* port 1, 4 bytes 12 */

            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps  rA1, rB1             /* port 0, 3 bytes 16 */
            addps  rB1, rC11            /* port 1, 4 bytes 04 */
            movaps 4192(pB), rB1         /* port 2, 7 bytes 11 */
            pshufd $0xAA, rB2, rB0      /* port 5, 5 bytes 16 */

            movaps rA0, rm0             /* port 5, 3 bytes 03 */
            mulps rB0, rm0              /* port 0, 3 bytes 06 */
            addps rm0, rC02             /* port 1, 4 bytes 10 */

            .byte 0x3e                  /* no prt, 1 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            .byte 0x3e                  /* no prt, 1 bytes 13 */
            mulps rA1, rB0              /* port 0, 3 bytes 16 */
            addps rB0, rC12             /* port 1, 4 bytes 04 */
            pshufd $0xFF, rB2, rm0      /* port 5, 5 bytes 09 */


            mulps rm0, rA0              /* port 0, 3 bytes 12 */
            addps rA0, rC03             /* port 1, 4 bytes 16 */
            pshufd $0x00, rB1, rB0      /* port 5, 5 bytes 05 */

            mulps rm0, rA1              /* port 0, 3 bytes 08 */
            movaps rB1, rB2             /* port 5, 3 bytes 11 */
            .byte 0x3e                  /* no prt, 1 bytes 12 */
            addps rA1, rC13             /* port 1, 4 bytes 16 */
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
         add incPF, pfA
         add incPF, pfB
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
#if 0
.global findSize
findSize:
mov $SS1-SS0, %rax
ret
SS0:
SS1:
#endif
