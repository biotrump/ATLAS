#include "atlas_asm.h"

#define rB0     %xmm0
#define rC00    %xmm1
#define rC01    %xmm2
#define rC02    %xmm3
#define rC03    %xmm4
#define rC04    %xmm5
#define rC05    %xmm6
#define rC06    %xmm7
#define rC07    %xmm8
#define rC08    %xmm9
#define rC09    %xmm10
#define rC10    %xmm11
#define rC11    %xmm12
#define rC12    %xmm13
#define rC13    %xmm14
#define rC14    %xmm15
/*
 * Prioritize original registers for inner-loop operations, but inc regs
 * can be anything w/o changing opcode size, so use new regs for those
 */
#define KK      %rdx  /* API reg */
#define pA      %rcx  /* API reg */
#define pB      %rax  /* comes in as r9 */
#define r256    %r9   /* set after mov r9 to pC () */
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
   sub $-128, pA
   sub $-128, pB
   sub $-128, pC
   sub $-128, pfB
   sub $-128, pfA
   mov $256, r256
   mov pA, pA0
   mov pB, pB0
   mov $15*1*8, incPF
/*
 * incAm = 30*sizeof*K = 30*8*K = 16*15*K
 */
   lea (KK, KK,4), incAm        /* incAm = 5*K */
   lea (incAm, incAm,2), incAm  /* incAM = 15*K */
   shl $4, incAm                /* incAm = 16*15*K */

   ALIGN16
   MLOOP:
      NLOOP:
/*
 * pB: 4 (add const) + 4 = 8 bytes
 * FMAC: 6*#
 * pA: 3 (add reg)
 * sz = 11 + 6*14 = 11 + 84 = 95  -> 6.79 bytes/FMAC
 * sz15 = 11 + 6*15 = 101  -> 6.733 bytes/FMAC
 *
 */
/*
 *       First iteration peeled to handle init of rC
 */
         vmovddup -128(pB), rB0
         vmulpd -128(pA), rB0, rC00
         vmulpd -112(pA), rB0, rC01
         vmulpd -96(pA), rB0, rC02
         vmulpd -80(pA), rB0, rC03
         vmulpd -64(pA), rB0, rC04
         vmulpd -48(pA), rB0, rC05
         vmulpd -32(pA), rB0, rC06
         vmulpd -16(pA), rB0, rC07
         vmulpd (pA), rB0, rC08
         vmulpd 16(pA), rB0, rC09
         vmulpd 32(pA), rB0, rC10
         vmulpd 48(pA), rB0, rC11
         vmulpd 64(pA), rB0, rC12
         vmulpd 80(pA), rB0, rC13
         vmulpd 96(pA), rB0, rC14

/*
 *       ==========================
 *       Completely unrolled K-loop
 *       ==========================
 */
         #if KB > 1
            vmovddup -120(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 2
            vmovddup -112(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 3
            vmovddup -104(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 4
            vmovddup -96(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 5
            vmovddup -88(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 6
            vmovddup -80(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 7
            vmovddup -72(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 8
            vmovddup -64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 9
            vmovddup -56(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 10
            vmovddup -48(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 11
            vmovddup -40(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 12
            vmovddup -32(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 13
            vmovddup -24(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 14
            vmovddup -16(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 15
            vmovddup -8(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 16
            vmovddup 0(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 17
            vmovddup 8(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 18
            vmovddup 16(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 19
            vmovddup 24(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 20
            vmovddup 32(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 21
            vmovddup 40(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 22
            vmovddup 48(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 23
            vmovddup 56(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 24
            vmovddup 64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 25
            vmovddup 72(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 26
            vmovddup 80(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 27
            vmovddup 88(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 28
            vmovddup 96(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 29
            vmovddup 104(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 30
            vmovddup 112(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 31
            vmovddup 120(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            add r256, pB
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 32
            vmovddup -128(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 33
            vmovddup -120(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 34
            vmovddup -112(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 35
            vmovddup -104(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 36
            vmovddup -96(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 37
            vmovddup -88(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 38
            vmovddup -80(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 39
            vmovddup -72(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 40
            vmovddup -64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 41
            vmovddup -56(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 42
            vmovddup -48(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 43
            vmovddup -40(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 44
            vmovddup -32(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 45
            vmovddup -24(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 46
            vmovddup -16(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 47
            vmovddup -8(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 48
            vmovddup 0(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 49
            vmovddup 8(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 50
            vmovddup 16(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 51
            vmovddup 24(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 52
            vmovddup 32(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 53
            vmovddup 40(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 54
            vmovddup 48(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 55
            vmovddup 56(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 56
            vmovddup 64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 57
            vmovddup 72(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 58
            vmovddup 80(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 59
            vmovddup 88(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 60
            vmovddup 96(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 61
            vmovddup 104(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 62
            vmovddup 112(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 63
            vmovddup 120(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            add r256, pB
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 64
            vmovddup -128(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 65
            vmovddup -120(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 66
            vmovddup -112(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 67
            vmovddup -104(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 68
            vmovddup -96(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 69
            vmovddup -88(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 70
            vmovddup -80(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 71
            vmovddup -72(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 72
            vmovddup -64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 73
            vmovddup -56(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 74
            vmovddup -48(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 75
            vmovddup -40(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 76
            vmovddup -32(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 77
            vmovddup -24(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 78
            vmovddup -16(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 79
            vmovddup -8(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 80
            vmovddup 0(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 81
            vmovddup 8(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 82
            vmovddup 16(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 83
            vmovddup 24(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 84
            vmovddup 32(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 85
            vmovddup 40(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 86
            vmovddup 48(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 87
            vmovddup 56(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 88
            vmovddup 64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 89
            vmovddup 72(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 90
            vmovddup 80(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 91
            vmovddup 88(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 92
            vmovddup 96(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 93
            vmovddup 104(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 94
            vmovddup 112(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 95
            vmovddup 120(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            add r256, pB
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 96
            vmovddup -128(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 97
            vmovddup -120(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 98
            vmovddup -112(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 99
            vmovddup -104(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 100
            vmovddup -96(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 101
            vmovddup -88(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 102
            vmovddup -80(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 103
            vmovddup -72(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 104
            vmovddup -64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 105
            vmovddup -56(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 106
            vmovddup -48(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 107
            vmovddup -40(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 108
            vmovddup -32(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 109
            vmovddup -24(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 110
            vmovddup -16(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 111
            vmovddup -8(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 112
            vmovddup 0(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 113
            vmovddup 8(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 114
            vmovddup 16(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 115
            vmovddup 24(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 116
            vmovddup 32(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 117
            vmovddup 40(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 118
            vmovddup 48(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 119
            vmovddup 56(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 120
            vmovddup 64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 121
            vmovddup 72(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 122
            vmovddup 80(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 123
            vmovddup 88(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 124
            vmovddup 96(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 125
            vmovddup 104(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 126
            vmovddup 112(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 127
            vmovddup 120(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            add r256, pB
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 128
            vmovddup -128(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 129
            vmovddup -120(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 130
            vmovddup -112(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 131
            vmovddup -104(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 132
            vmovddup -96(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 133
            vmovddup -88(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 134
            vmovddup -80(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 135
            vmovddup -72(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 136
            vmovddup -64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 137
            vmovddup -56(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 138
            vmovddup -48(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 139
            vmovddup -40(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 140
            vmovddup -32(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 141
            vmovddup -24(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 142
            vmovddup -16(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 143
            vmovddup -8(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 144
            vmovddup 0(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 145
            vmovddup 8(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 146
            vmovddup 16(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 147
            vmovddup 24(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 148
            vmovddup 32(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 149
            vmovddup 40(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 150
            vmovddup 48(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 151
            vmovddup 56(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 152
            vmovddup 64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 153
            vmovddup 72(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 154
            vmovddup 80(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 155
            vmovddup 88(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 156
            vmovddup 96(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 157
            vmovddup 104(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 158
            vmovddup 112(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 159
            vmovddup 120(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            add r256, pB
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 160
            vmovddup -128(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 161
            vmovddup -120(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 162
            vmovddup -112(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 163
            vmovddup -104(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 164
            vmovddup -96(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 165
            vmovddup -88(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 166
            vmovddup -80(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 167
            vmovddup -72(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 168
            vmovddup -64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 169
            vmovddup -56(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 170
            vmovddup -48(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 171
            vmovddup -40(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 172
            vmovddup -32(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 173
            vmovddup -24(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 174
            vmovddup -16(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 175
            vmovddup -8(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 176
            vmovddup 0(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 177
            vmovddup 8(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 178
            vmovddup 16(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 179
            vmovddup 24(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 180
            vmovddup 32(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 181
            vmovddup 40(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 182
            vmovddup 48(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 183
            vmovddup 56(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 184
            vmovddup 64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 185
            vmovddup 72(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 186
            vmovddup 80(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 187
            vmovddup 88(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 188
            vmovddup 96(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 189
            vmovddup 104(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 190
            vmovddup 112(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 191
            vmovddup 120(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            add r256, pB
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 192
            vmovddup -128(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 193
            vmovddup -120(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 194
            vmovddup -112(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 195
            vmovddup -104(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 196
            vmovddup -96(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 197
            vmovddup -88(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 198
            vmovddup -80(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 199
            vmovddup -72(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 200
            vmovddup -64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 201
            vmovddup -56(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 202
            vmovddup -48(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 203
            vmovddup -40(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 204
            vmovddup -32(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 205
            vmovddup -24(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 206
            vmovddup -16(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 207
            vmovddup -8(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 208
            vmovddup 0(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 209
            vmovddup 8(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 210
            vmovddup 16(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 211
            vmovddup 24(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 212
            vmovddup 32(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 213
            vmovddup 40(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 214
            vmovddup 48(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 215
            vmovddup 56(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 216
            vmovddup 64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 217
            vmovddup 72(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 218
            vmovddup 80(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 219
            vmovddup 88(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 220
            vmovddup 96(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 221
            vmovddup 104(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 222
            vmovddup 112(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 223
            vmovddup 120(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            add r256, pB
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 224
            vmovddup -128(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 225
            vmovddup -120(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 226
            vmovddup -112(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 227
            vmovddup -104(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 228
            vmovddup -96(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 229
            vmovddup -88(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 230
            vmovddup -80(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 231
            vmovddup -72(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 232
            vmovddup -64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 233
            vmovddup -56(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 234
            vmovddup -48(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 235
            vmovddup -40(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 236
            vmovddup -32(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 237
            vmovddup -24(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 238
            vmovddup -16(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 239
            vmovddup -8(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif
         #if KB > 240
            vmovddup 0(pB), rB0
            FMAC -128(pA), rB0, rC00
            FMAC -112(pA), rB0, rC01
            FMAC -96(pA), rB0, rC02
            FMAC -80(pA), rB0, rC03
            FMAC -64(pA), rB0, rC04
            FMAC -48(pA), rB0, rC05
            FMAC -32(pA), rB0, rC06
            FMAC -16(pA), rB0, rC07
            FMAC 0(pA), rB0, rC08
            FMAC 16(pA), rB0, rC09
            FMAC 32(pA), rB0, rC10
            FMAC 48(pA), rB0, rC11
            FMAC 64(pA), rB0, rC12
            FMAC 80(pA), rB0, rC13
            FMAC 96(pA), rB0, rC14
         #endif
         #if KB > 241
            vmovddup 8(pB), rB0
            FMAC 112(pA), rB0, rC00
            add r256, pA
            FMAC -128(pA), rB0, rC01
            FMAC -112(pA), rB0, rC02
            FMAC -96(pA), rB0, rC03
            FMAC -80(pA), rB0, rC04
            FMAC -64(pA), rB0, rC05
            FMAC -48(pA), rB0, rC06
            FMAC -32(pA), rB0, rC07
            FMAC -16(pA), rB0, rC08
            FMAC 0(pA), rB0, rC09
            FMAC 16(pA), rB0, rC10
            FMAC 32(pA), rB0, rC11
            FMAC 48(pA), rB0, rC12
            FMAC 64(pA), rB0, rC13
            FMAC 80(pA), rB0, rC14
         #endif
         #if KB > 242
            vmovddup 16(pB), rB0
            FMAC 96(pA), rB0, rC00
            FMAC 112(pA), rB0, rC01
            add r256, pA
            FMAC -128(pA), rB0, rC02
            FMAC -112(pA), rB0, rC03
            FMAC -96(pA), rB0, rC04
            FMAC -80(pA), rB0, rC05
            FMAC -64(pA), rB0, rC06
            FMAC -48(pA), rB0, rC07
            FMAC -32(pA), rB0, rC08
            FMAC -16(pA), rB0, rC09
            FMAC 0(pA), rB0, rC10
            FMAC 16(pA), rB0, rC11
            FMAC 32(pA), rB0, rC12
            FMAC 48(pA), rB0, rC13
            FMAC 64(pA), rB0, rC14
         #endif
         #if KB > 243
            vmovddup 24(pB), rB0
            FMAC 80(pA), rB0, rC00
            FMAC 96(pA), rB0, rC01
            FMAC 112(pA), rB0, rC02
            add r256, pA
            FMAC -128(pA), rB0, rC03
            FMAC -112(pA), rB0, rC04
            FMAC -96(pA), rB0, rC05
            FMAC -80(pA), rB0, rC06
            FMAC -64(pA), rB0, rC07
            FMAC -48(pA), rB0, rC08
            FMAC -32(pA), rB0, rC09
            FMAC -16(pA), rB0, rC10
            FMAC 0(pA), rB0, rC11
            FMAC 16(pA), rB0, rC12
            FMAC 32(pA), rB0, rC13
            FMAC 48(pA), rB0, rC14
         #endif
         #if KB > 244
            vmovddup 32(pB), rB0
            FMAC 64(pA), rB0, rC00
            FMAC 80(pA), rB0, rC01
            FMAC 96(pA), rB0, rC02
            FMAC 112(pA), rB0, rC03
            add r256, pA
            FMAC -128(pA), rB0, rC04
            FMAC -112(pA), rB0, rC05
            FMAC -96(pA), rB0, rC06
            FMAC -80(pA), rB0, rC07
            FMAC -64(pA), rB0, rC08
            FMAC -48(pA), rB0, rC09
            FMAC -32(pA), rB0, rC10
            FMAC -16(pA), rB0, rC11
            FMAC 0(pA), rB0, rC12
            FMAC 16(pA), rB0, rC13
            FMAC 32(pA), rB0, rC14
         #endif
         #if KB > 245
            vmovddup 40(pB), rB0
            FMAC 48(pA), rB0, rC00
            FMAC 64(pA), rB0, rC01
            FMAC 80(pA), rB0, rC02
            FMAC 96(pA), rB0, rC03
            FMAC 112(pA), rB0, rC04
            add r256, pA
            FMAC -128(pA), rB0, rC05
            FMAC -112(pA), rB0, rC06
            FMAC -96(pA), rB0, rC07
            FMAC -80(pA), rB0, rC08
            FMAC -64(pA), rB0, rC09
            FMAC -48(pA), rB0, rC10
            FMAC -32(pA), rB0, rC11
            FMAC -16(pA), rB0, rC12
            FMAC 0(pA), rB0, rC13
            FMAC 16(pA), rB0, rC14
         #endif
         #if KB > 246
            vmovddup 48(pB), rB0
            FMAC 32(pA), rB0, rC00
            FMAC 48(pA), rB0, rC01
            FMAC 64(pA), rB0, rC02
            FMAC 80(pA), rB0, rC03
            FMAC 96(pA), rB0, rC04
            FMAC 112(pA), rB0, rC05
            add r256, pA
            FMAC -128(pA), rB0, rC06
            FMAC -112(pA), rB0, rC07
            FMAC -96(pA), rB0, rC08
            FMAC -80(pA), rB0, rC09
            FMAC -64(pA), rB0, rC10
            FMAC -48(pA), rB0, rC11
            FMAC -32(pA), rB0, rC12
            FMAC -16(pA), rB0, rC13
            FMAC 0(pA), rB0, rC14
         #endif
         #if KB > 247
            vmovddup 56(pB), rB0
            FMAC 16(pA), rB0, rC00
            FMAC 32(pA), rB0, rC01
            FMAC 48(pA), rB0, rC02
            FMAC 64(pA), rB0, rC03
            FMAC 80(pA), rB0, rC04
            FMAC 96(pA), rB0, rC05
            FMAC 112(pA), rB0, rC06
            add r256, pA
            FMAC -128(pA), rB0, rC07
            FMAC -112(pA), rB0, rC08
            FMAC -96(pA), rB0, rC09
            FMAC -80(pA), rB0, rC10
            FMAC -64(pA), rB0, rC11
            FMAC -48(pA), rB0, rC12
            FMAC -32(pA), rB0, rC13
            FMAC -16(pA), rB0, rC14
         #endif
         #if KB > 248
            vmovddup 64(pB), rB0
            FMAC 0(pA), rB0, rC00
            FMAC 16(pA), rB0, rC01
            FMAC 32(pA), rB0, rC02
            FMAC 48(pA), rB0, rC03
            FMAC 64(pA), rB0, rC04
            FMAC 80(pA), rB0, rC05
            FMAC 96(pA), rB0, rC06
            FMAC 112(pA), rB0, rC07
            add r256, pA
            FMAC -128(pA), rB0, rC08
            FMAC -112(pA), rB0, rC09
            FMAC -96(pA), rB0, rC10
            FMAC -80(pA), rB0, rC11
            FMAC -64(pA), rB0, rC12
            FMAC -48(pA), rB0, rC13
            FMAC -32(pA), rB0, rC14
         #endif
         #if KB > 249
            vmovddup 72(pB), rB0
            FMAC -16(pA), rB0, rC00
            FMAC 0(pA), rB0, rC01
            FMAC 16(pA), rB0, rC02
            FMAC 32(pA), rB0, rC03
            FMAC 48(pA), rB0, rC04
            FMAC 64(pA), rB0, rC05
            FMAC 80(pA), rB0, rC06
            FMAC 96(pA), rB0, rC07
            FMAC 112(pA), rB0, rC08
            add r256, pA
            FMAC -128(pA), rB0, rC09
            FMAC -112(pA), rB0, rC10
            FMAC -96(pA), rB0, rC11
            FMAC -80(pA), rB0, rC12
            FMAC -64(pA), rB0, rC13
            FMAC -48(pA), rB0, rC14
         #endif
         #if KB > 250
            vmovddup 80(pB), rB0
            FMAC -32(pA), rB0, rC00
            FMAC -16(pA), rB0, rC01
            FMAC 0(pA), rB0, rC02
            FMAC 16(pA), rB0, rC03
            FMAC 32(pA), rB0, rC04
            FMAC 48(pA), rB0, rC05
            FMAC 64(pA), rB0, rC06
            FMAC 80(pA), rB0, rC07
            FMAC 96(pA), rB0, rC08
            FMAC 112(pA), rB0, rC09
            add r256, pA
            FMAC -128(pA), rB0, rC10
            FMAC -112(pA), rB0, rC11
            FMAC -96(pA), rB0, rC12
            FMAC -80(pA), rB0, rC13
            FMAC -64(pA), rB0, rC14
         #endif
         #if KB > 251
            vmovddup 88(pB), rB0
            FMAC -48(pA), rB0, rC00
            FMAC -32(pA), rB0, rC01
            FMAC -16(pA), rB0, rC02
            FMAC 0(pA), rB0, rC03
            FMAC 16(pA), rB0, rC04
            FMAC 32(pA), rB0, rC05
            FMAC 48(pA), rB0, rC06
            FMAC 64(pA), rB0, rC07
            FMAC 80(pA), rB0, rC08
            FMAC 96(pA), rB0, rC09
            FMAC 112(pA), rB0, rC10
            add r256, pA
            FMAC -128(pA), rB0, rC11
            FMAC -112(pA), rB0, rC12
            FMAC -96(pA), rB0, rC13
            FMAC -80(pA), rB0, rC14
         #endif
         #if KB > 252
            vmovddup 96(pB), rB0
            FMAC -64(pA), rB0, rC00
            FMAC -48(pA), rB0, rC01
            FMAC -32(pA), rB0, rC02
            FMAC -16(pA), rB0, rC03
            FMAC 0(pA), rB0, rC04
            FMAC 16(pA), rB0, rC05
            FMAC 32(pA), rB0, rC06
            FMAC 48(pA), rB0, rC07
            FMAC 64(pA), rB0, rC08
            FMAC 80(pA), rB0, rC09
            FMAC 96(pA), rB0, rC10
            FMAC 112(pA), rB0, rC11
            add r256, pA
            FMAC -128(pA), rB0, rC12
            FMAC -112(pA), rB0, rC13
            FMAC -96(pA), rB0, rC14
         #endif
         #if KB > 253
            vmovddup 104(pB), rB0
            FMAC -80(pA), rB0, rC00
            FMAC -64(pA), rB0, rC01
            FMAC -48(pA), rB0, rC02
            FMAC -32(pA), rB0, rC03
            FMAC -16(pA), rB0, rC04
            FMAC 0(pA), rB0, rC05
            FMAC 16(pA), rB0, rC06
            FMAC 32(pA), rB0, rC07
            FMAC 48(pA), rB0, rC08
            FMAC 64(pA), rB0, rC09
            FMAC 80(pA), rB0, rC10
            FMAC 96(pA), rB0, rC11
            FMAC 112(pA), rB0, rC12
            add r256, pA
            FMAC -128(pA), rB0, rC13
            FMAC -112(pA), rB0, rC14
         #endif
         #if KB > 254
            vmovddup 112(pB), rB0
            FMAC -96(pA), rB0, rC00
            FMAC -80(pA), rB0, rC01
            FMAC -64(pA), rB0, rC02
            FMAC -48(pA), rB0, rC03
            FMAC -32(pA), rB0, rC04
            FMAC -16(pA), rB0, rC05
            FMAC 0(pA), rB0, rC06
            FMAC 16(pA), rB0, rC07
            FMAC 32(pA), rB0, rC08
            FMAC 48(pA), rB0, rC09
            FMAC 64(pA), rB0, rC10
            FMAC 80(pA), rB0, rC11
            FMAC 96(pA), rB0, rC12
            FMAC 112(pA), rB0, rC13
            add r256, pA
            FMAC -128(pA), rB0, rC14
         #endif
         #if KB > 255
            vmovddup 120(pB), rB0
            FMAC -112(pA), rB0, rC00
            FMAC -96(pA), rB0, rC01
            FMAC -80(pA), rB0, rC02
            FMAC -64(pA), rB0, rC03
            add r256, pB
            FMAC -48(pA), rB0, rC04
            FMAC -32(pA), rB0, rC05
            FMAC -16(pA), rB0, rC06
            FMAC 0(pA), rB0, rC07
            FMAC 16(pA), rB0, rC08
            FMAC 32(pA), rB0, rC09
            FMAC 48(pA), rB0, rC10
            FMAC 64(pA), rB0, rC11
            FMAC 80(pA), rB0, rC12
            FMAC 96(pA), rB0, rC13
            FMAC 112(pA), rB0, rC14
            add r256, pA
         #endif

         #ifndef BETA0
            VCOP -128(pC), rC00
         #endif
         movapd rC00, -128(pC)
         #ifndef BETA0
            VCOP -112(pC), rC01
         #endif
         movapd rC01, -112(pC)
         #ifndef BETA0
            VCOP -96(pC), rC02
         #endif
         movapd rC02, -96(pC)
         #ifndef BETA0
            VCOP -80(pC), rC03
         #endif
         movapd rC03, -80(pC)
         #ifndef BETA0
            VCOP -64(pC), rC04
         #endif
         movapd rC04, -64(pC)
         #ifndef BETA0
            VCOP -48(pC), rC05
         #endif
         movapd rC05, -48(pC)
         #ifndef BETA0
            VCOP -32(pC), rC06
         #endif
         movapd rC06, -32(pC)
         #ifndef BETA0
            VCOP -16(pC), rC07
         #endif
         movapd rC07, -16(pC)
         #ifndef BETA0
            VCOP (pC), rC08
         #endif
         movapd rC08, (pC)
         #ifndef BETA0
            VCOP 16(pC), rC09
         #endif
         movapd rC09, 16(pC)
         #ifndef BETA0
            VCOP 32(pC), rC10
         #endif
         movapd rC10, 32(pC)
         #ifndef BETA0
            VCOP 48(pC), rC11
         #endif
         movapd rC11, 48(pC)
         #ifndef BETA0
            VCOP 64(pC), rC12
         #endif
         movapd rC12, 64(pC)
         #ifndef BETA0
            VCOP 80(pC), rC13
         #endif
         movapd rC13, 80(pC)
         #ifndef BETA0
            VCOP 96(pC), rC14
         #endif
         movapd rC14, 96(pC)
         add $240, pC
         mov pA0, pA
         add $KB*8-(((KB*8)/256)*256), pB
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
