#include "atlas_asm.h"

#define rB0     %xmm0
#define rB1     %xmm1
#define rB2     %xmm2
#define rA0     %xmm3
#define rA1     %xmm4
#define rA2     %xmm5
#define rC00    %xmm6
#define rC10    %xmm7
#define rC20    %xmm8
#define rC01    %xmm9
#define rC11    %xmm10
#define rC21    %xmm11
#define rC02    %xmm12
#define rC12    %xmm13
#define rC22    %xmm14
/*
 * Prioritize original registers for inner-loop operations, but inc regs
 * can be anything w/o changing opcode size, so use new regs for those
 */
#define KK      %rdx  /* API reg */
#define pA      %rcx  /* API reg */
#define pB      %rax  /* comes in as r9 */
#define r24     %r9   /* set after mov r9 to pC () */
/*
 * Then N-loop variables much less important, so use any orig regs left
 */
#define pA0     %r8   /* set after mov r8 to pB (rax) */
#define pC      %rsi  /* set after mov rsi to nnu () */
#define nnu     %r10  /* comes in as rsi */
#define pfA     %rbx
#define pfB     %rbp
#define incPF   %r12
#define KK0     %rdi
/*
 * We could give a rat's ass about what registers used in outer (M-) loop
 */
#define nmu     %r11  /* comes in as rdi */
#define incAm   %r13
#define nnu0    %r14
#define pB0     %r15
/*
                    rdi      rsi    rdx        rcx         r8        r9
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                  8(%rsp)    16(%rsp)     24(%rsp)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */
#define PFBDIST 128
#define PFADIST 128
#define prefA(m_) prefetcht0 m_
#define prefB(m_) prefetcht0 m_
#define prefC(m_) prefetcht0 m_
#define FMAC vfmadd231pd   /* FMAC m256/r256, rs1, rd */
#if defined(BETAN) || defined(BETAn)
   #define BETAN1
#endif
#ifdef BETAN1
   #define VCOP subpd
#else
   #define VCOP addpd
#endif
#define movapd movaps
.text
.global ATL_asmdecor(ATL_USERMM)
ALIGN16
ATL_asmdecor(ATL_USERMM):
/*
 * Save callee-saved iregs
 */
   movq %rbp, -8(%rsp)
   movq %rbx, -16(%rsp)
   movq %r12, -24(%rsp)
   movq %r13, -32(%rsp)
   movq %r14, -40(%rsp)
   movq %r15, -48(%rsp)
/*
 * Load paramaters
 */
   mov %rdi, nmu
   mov %rsi, nnu
   mov %r8, pB
   mov %r9, pC
   mov nnu, nnu0
   movq 8(%rsp), pfB      /* pfB = pAn */
   movq 16(%rsp), pfA     /* pf = pBn */
   cmp pfB, pB
   CMOVE pfA, pfB
   mov KK, KK0
   sub $-128, pC
   sub $-128, pfB
   sub $-128, pfA
   mov $24, r24
   mov pA, pA0
   mov pB, pB0
   mov $6*3*8, incPF
/*
 * incAm = 6*sizeof*K = 6*8*K = 3*2*8*K = 3*K*16
 */
   lea (KK, KK,2), incAm   /* incAm = 3*K */
   shl $4, incAm           /* incAm = 32*3*K */
   vxorps rA1, rA1, rA1
   vxorps rA2, rA2, rA2

   ALIGN16
   MLOOP:
      NLOOP:
/*
 *       First iteration peeled to handle init of rC
 */
         vmovddup (pB), rB0     /* 4  04 */
         vmovaps (pA), rA0      /* 4  08 */
         vmulpd rA0, rB0, rC00  /* 4  12 */
         add $48, pA            /* 4  16 */

         movaps -32(pA), rA1    /* 4  04 */
         vmulpd rA1, rB0, rC10  /* 4  08 */
         movaps -16(pA), rA2    /* 4  12 */
         vmulpd rA2, rB0, rC20  /* 4  16 */

         vmovddup 8(pB), rB1    /* 5  05 */
         vmulpd rA0, rB1, rC01  /* 4  09 */
            prefC(-128(pC))     /* 4  13 */
         add r24, pB            /* 3  16 */

         vmulpd rA1, rB1, rC11  /* 4  04 */
            .byte 0x3e
            prefB(-128(pfB))    /* 4  08 */
         vmulpd rA2, rB1, rC21  /* 4  12 */
            add incPF, pfB      /* 3  15 */

         vmovddup -8(pB), rB2   /* 5  05 */
         vmulpd rA0, rB2, rC02  /* 4  09 */
            prefC((pC))         /* 3  12 */
         vmulpd rA1, rB2, rC12  /* 4  16 */

         vmulpd rA2, rB2, rC22  /* 4  04 */
            prefB(-16(pfB))     /* 4  08 */
            .byte 0x3e
            prefA(-128(pfA))    /* 4  12 */
            add incPF, pfA      /* 3  15 */
/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 1
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 2
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 3
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 4
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 5
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 6
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 7
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 8
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 9
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 10
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 11
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 12
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 13
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 14
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 15
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 16
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 17
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 18
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 19
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 20
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 21
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 22
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 23
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 24
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 25
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 26
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 27
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 28
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 29
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 30
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 31
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 32
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 33
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 34
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 35
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 36
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 37
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 38
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 39
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 40
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 41
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 42
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 43
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 44
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 45
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 46
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 47
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 48
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 49
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 50
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 51
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 52
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 53
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 54
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 55
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 56
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 57
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 58
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 59
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 60
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 61
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 62
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 63
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 64
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 65
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 66
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 67
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 68
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 69
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 70
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 71
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 72
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 73
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 74
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 75
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 76
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 77
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 78
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 79
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 80
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 81
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 82
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 83
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 84
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 85
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 86
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 87
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 88
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 89
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 90
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 91
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 92
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 93
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 94
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 95
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 96
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 97
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 98
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 99
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 100
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 101
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 102
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 103
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 104
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 105
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 106
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 107
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 108
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 109
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 110
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 111
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 112
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 113
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 114
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 115
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 116
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 117
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 118
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 119
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 120
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 121
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 122
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 123
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 124
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 125
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 126
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 127
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 128
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 129
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 130
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 131
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 132
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 133
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 134
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 135
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 136
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 137
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 138
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 139
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 140
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 141
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 142
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 143
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 144
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 145
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 146
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 147
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 148
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 149
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 150
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 151
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 152
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 153
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 154
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 155
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 156
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 157
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 158
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 159
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 160
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 161
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 162
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 163
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 164
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 165
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 166
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 167
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 168
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 169
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 170
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 171
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 172
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 173
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 174
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 175
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 176
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 177
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 178
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 179
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 180
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 181
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 182
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 183
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 184
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 185
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 186
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 187
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 188
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 189
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 190
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 191
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 192
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 193
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 194
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 195
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 196
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 197
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 198
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 199
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 200
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 201
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 202
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 203
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 204
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 205
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 206
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 207
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 208
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 209
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 210
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 211
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 212
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 213
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 214
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 215
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 216
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 217
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 218
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 219
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 220
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 221
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 222
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 223
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 224
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 225
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 226
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 227
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 228
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 229
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 230
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 231
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 232
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 233
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 234
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 235
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 236
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 237
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 238
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 239
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 240
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 241
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 242
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 243
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 244
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 245
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 246
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 247
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 248
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 249
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 250
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 251
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 252
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 253
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 254
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif
         #if KB > 255
            movaps (pA), rA0               /* 3 03 */
            movddup (pB), rB0              /* 4 04 */
            FMAC rA0, rB0, rC00            /* 5 12 */
            add $48, pA                    /* 4 16 */

            movaps -32(pA), rA1            /* 4 04 */
            FMAC rA1, rB0, rC10            /* 5 09 */
            movaps -16(pA), rA2            /* 4 13 */
            add r24, pB                    /* 3 16 */

            FMAC rA2, rB0, rC20            /* 5 05 */
            .byte 0x3e
            movddup -16(pB), rB1           /* 5 10 */
            FMAC rA0, rB1, rC01            /* 5 15 */

            FMAC rA1, rB1, rC11            /* 5 05 */
            .byte 0x3e
            movddup -8(pB), rB2            /* 5 10 */
            FMAC rA2, rB1, rC21            /* 5 15 */

            FMAC rA0, rB2, rC02            /* 5 05 */
            FMAC rA1, rB2, rC12            /* 5 05 */
            .byte 0x3e
            FMAC rA2, rB2, rC22            /* 5 05 */
         #endif

         #ifndef BETA0
            VCOP -128(pC), rC00
         #endif
         movapd rC00, -128(pC)
         #ifndef BETA0
            VCOP -112(pC), rC10
         #endif
         movapd rC10, -112(pC)
         #ifndef BETA0
            VCOP -96(pC), rC20
         #endif
         movapd rC20, -96(pC)
         #ifndef BETA0
            VCOP -80(pC), rC01
         #endif
         movapd rC01, -80(pC)
         #ifndef BETA0
            VCOP -64(pC), rC11
         #endif
         movapd rC11, -64(pC)
         #ifndef BETA0
            VCOP -48(pC), rC21
         #endif
         movapd rC21, -48(pC)
         #ifndef BETA0
            VCOP -32(pC), rC02
         #endif
         movapd rC02, -32(pC)
         #ifndef BETA0
            VCOP -16(pC), rC12
         #endif
         movapd rC12, -16(pC)
         #ifndef BETA0
            VCOP (pC), rC22
         #endif
         movapd rC22, (pC)
         mov pA0, pA
         add $144, pC
         sub $1, nnu
      jnz NLOOP
      mov nnu0, nnu
      add incAm, pA0
      mov pA0, pA
      mov pB0, pB
      sub $1, nmu
   jnz MLOOP
 DONE:
   movq -8(%rsp), %rbp
   movq -16(%rsp), %rbx
   movq -24(%rsp), %r12
   movq -32(%rsp), %r13
   movq -40(%rsp), %r14
   movq -48(%rsp), %r15
   ret
#if 0
.global findSize
findSize:
mov $SS1-SS0, %rax
ret
SS0:
SS1:
#endif
