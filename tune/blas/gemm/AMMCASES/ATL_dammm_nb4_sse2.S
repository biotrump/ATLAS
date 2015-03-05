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

#define pA %rcx
#define pB %r8
#define pC %r9
#define m0      %xmm0
#define rA00    %xmm1
#define rA10    %xmm2
#define rA01    %xmm3
#define rA11    %xmm4
#define rA02    %xmm5
#define rA12    %xmm6
#define rA03    %xmm7
#define rA13    %xmm8
#define rB0     %xmm9
#define rB1     %xmm10
#define rC00    %xmm11
#define rC10    %xmm12
#define rC01    %xmm13
#define rC11    %xmm14
/*
                    rdi      rsi    rdx        rcx         r8        r9
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                  8(%rsp)    16(%rsp)     24(%rsp)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */
.text
ALIGN16
.global ATL_asmdecor(ATL_USERMM)
ATL_asmdecor(ATL_USERMM):
   movapd (pA), rA00
   movapd 16(pA), rA10
   movapd 32(pA), rA01
   movapd 48(pA), rA11
   movapd 64(pA), rA02
   movapd 80(pA), rA12
   movapd 96(pA), rA03
   movapd 112(pA), rA13

/*
 * ======================================================
 * rC00 = rA00 * rB00 + rA01*rB10 + rA02*rB20 + rA03*rB30
 * rC10 = rA10 * rB00 + rA11*rB10 + rA12*rB20 + rA13*rB30
 * rC01 = rA00 * rB01 + rA01*rB11 + rA02*rB21 + rA03*rB31
 * rC11 = rA10 * rB01 + rA11*rB11 + rA12*rB21 + rA13*rB31
 * ======================================================
 */
/*
 * rC00 = rA00 * rB0
 * rC10 = rA10 * rB0
 * rC01 = rA00 * rB1
 * rC11 = rA10 * rB1
 */
   #if defined(BETA0) || defined(BETAN1)
      movddup (pB), rC10
      movapd rC10, rC00
      mulpd rA00, rC00
      movddup 32(pB), rC11
      mulpd rA10, rC10
      movapd rC11, rC01
      mulpd rA00, rC01
      mulpd rA10, rC11
   #else
      movddup (pB), rB0
      movddup 32(pB), rB1
      movapd rB0, m0
      mulpd rA00, m0
      movapd (pC), rC00
      addpd m0, rC00
      mulpd rA10, rB0
      movapd 16(pC), rC10
      addpd rB0, rC10
      movapd rB1, m0
      mulpd rA00, m0
      movapd 32(pC), rC01
      addpd m0, rC01
      mulpd rA10, rB1
      movapd 48(pC), rC11
      addpd rB1, rC11
   #endif
/*
 * rC00 = rA01*rB10
 * rC10 = rA11*rB10
 * rC01 = rA01*rB11
 * rC11 = rA11*rB11
 */
   movddup 8(pB), rB0
   movapd rB0, m0
   mulpd rA01, m0
   addpd m0, rC00
   mulpd rA11, rB0
   addpd rB0, rC10
   movddup 40(pB), rB1
   movapd rB1, m0
   mulpd rA01, m0
   addpd m0, rC01
   mulpd rA11, rB1
   addpd rB1, rC11
/*
 * rC00 = rA02*rB20
 * rC10 = rA12*rB20
 * rC01 = rA02*rB21
 * rC11 = rA12*rB21
 */
   movddup 16(pB), rB0
   movapd rA02, m0
   mulpd rB0, m0
   addpd m0, rC00
   mulpd rA12, rB0
   addpd rB0, rC10
   movapd rA02, m0
   movddup 48(pB), rB1
   mulpd rB1, m0
   addpd m0, rC01
   mulpd rA12, rB1
   addpd rB1, rC11
/*
 * rC00 = rA03*rB30
 * rC10 = rA13*rB30
 * rC01 = rA03*rB31
 * rC11 = rA13*rB31
 */
   movddup 24(pB), rB0
   movapd rA03, m0
   mulpd rB0, m0
   addpd m0, rC00
   #ifdef BETAN1
      subpd (pC), rC00
   #endif
   movapd rC00, (pC)
   mulpd rA13, rB0
   addpd rB0, rC10
   #ifdef BETAN1
      subpd 16(pC), rC10
   #endif
   movapd rC10, 16(pC)
   movapd rA03, m0
   movddup 56(pB), rB1
   mulpd rB1, m0
   addpd m0, rC01
   #ifdef BETAN1
      subpd 32(pC), rC01
   #endif
   movapd rC01, 32(pC)
   mulpd rA13, rB1
   addpd rB1, rC11
   #ifdef BETAN1
      subpd 48(pC), rC11
   #endif
   movapd rC11, 48(pC)

/*
 * rC00 = rA00 * rB0
 * rC10 = rA10 * rB0
 * rC01 = rA00 * rB1
 * rC11 = rA10 * rB1
 */
   #if defined(BETA0) || defined(BETAN1)
      movddup 64(pB), rC10
      movapd rC10, rC00
      mulpd rA00, rC00
      movddup 96(pB), rC11
      mulpd rA10, rC10
      movapd rC11, rC01
      mulpd rA00, rC01
      mulpd rA10, rC11
   #else
      movddup 64(pB), rB0
      movddup 96(pB), rB1
      movapd rB0, m0
      mulpd rA00, m0
      movapd 64(pC), rC00
      addpd m0, rC00
      mulpd rA10, rB0
      movapd 80(pC), rC10
      addpd rB0, rC10
      movapd rB1, m0
      mulpd rA00, m0
      movapd 96(pC), rC01
      addpd m0, rC01
      mulpd rA10, rB1
      movapd 112(pC), rC11
      addpd rB1, rC11
   #endif
/*
 * rC00 = rA01*rB10
 * rC10 = rA11*rB10
 * rC01 = rA01*rB11
 * rC11 = rA11*rB11
 */
   movddup 72(pB), rB0
   movapd rB0, m0
   mulpd rA01, m0
   addpd m0, rC00
   mulpd rA11, rB0
   addpd rB0, rC10
   movddup 104(pB), rB1
   movapd rB1, m0
   mulpd rB1, rA01
   addpd rA01, rC01
   mulpd rA11, rB1
   addpd rB1, rC11
/*
 * rC00 = rA02*rB20
 * rC10 = rA12*rB20
 * rC01 = rA02*rB21
 * rC11 = rA12*rB21
 */
   movddup 80(pB), rB0
   movapd rA02, m0
   mulpd rB0, m0
   addpd m0, rC00
   mulpd rA12, rB0
   addpd rB0, rC10
   movddup 112(pB), rB1
   mulpd rB1, rA02
   addpd rA02, rC01
   mulpd rA12, rB1
   addpd rB1, rC11
/*
 * rC00 = rA03*rB30
 * rC10 = rA13*rB30
 * rC01 = rA03*rB31
 * rC11 = rA13*rB31
 */
   movddup 88(pB), rB0
   movapd rA03, m0
   mulpd rB0, m0
   addpd m0, rC00
   #ifdef BETAN1
      subpd 64(pC), rC00
   #endif
   movapd rC00, 64(pC)
   mulpd rA13, rB0
   addpd rB0, rC10
   #ifdef BETAN1
      subpd 80(pC), rC10
   #endif
   movapd rC10, 80(pC)
   movddup 120(pB), rB1
   mulpd rB1, rA03
   addpd rA03, rC01
   #ifdef BETAN1
      subpd 96(pC), rC01
   #endif
   movapd rC01, 96(pC)
   mulpd rA13, rB1
   addpd rB1, rC11
   #ifdef BETAN1
      subpd 112(pC), rC11
   #endif
   movapd rC11, 112(pC)
   ret
