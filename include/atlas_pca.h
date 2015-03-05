#ifndef ATLAS_PCA_H
   #define ATLAS_PCA_H

#include "atlas_misc.h"
#include "atlas_lapack.h"

/*
 * OpenMP provides horrible performance in general, but it is worse than serial
 * for PCA panel factorizations, so turn it off if the user has demanded OpenMP
 */
#ifdef ATL_OMP_THREADS
/*      #define ATL_USEPCA 1 */

/*
 * PowerPCs, POWERs and ARMs are weakly ordered, meaning that a given
 * processor's writes  may appear out-of-order to other processors,
 * which breaks PCA's syncs since PCA depends on in-order writes.
 * To fix, we must issue a memory barrier call before giving the go-ahead.
 * PowerPC: SYNC ensures that all prior stores complete before the next one.
 * POWER: DCS waits until all pending writes are written before preceeding
 * ARM: DMB (data mem barrier) - all prior mem accesses (in program order)
 *      complete before DMB returns
 *
 * Older x86's have a special mode where stores can become out-of-order, but
 * it was rarely enabled and does not seem to exist on modern hardware, so
 * we don't have to bother there.
 *
 * SPARCs do not change the order of stores.
 *
 * PowerPC and ARM syncs do not fix problem, so don't allow PCA on machines
 * with out-of-order write schemes.
 */
#elif defined(ATL_ARCH_PPCG4) || defined(ATL_ARCH_PPCG5)
   #ifdef __GNUC__
      #define ATL_membarrier __asm__ __volatile__ ("sync")
/*      #define ATL_USEPCA 1 */
   #endif
#elif defined(ATL_ARCH_POWER3) || defined(ATL_ARCH_POWER4) || \
      defined(ATL_ARCH_POWER5) || defined(ATL_ARCH_POWER6) || \
      defined(ATL_ARCH_POWER7)
   #ifdef __GNUC__
      #define ATL_membarrier __asm__ __volatile__ ("dcs")
/*      #define ATL_USEPCA 1 */
   #endif
/*
 * Unfortunately, none of the memory fence instructions seems to work
 * adequately on ARM
 */
#elif defined(ATL_ARCH_ARM64)
   #ifdef __GNUC__
      #define ATL_membarrier __asm__ __volatile__ ("dmb sy" : : : "memory")
/*      #define ATL_USEPCA 1 */
   #endif
#elif defined(ATL_ARCH_ARMv7)
   #ifdef __GNUC__
      #define ATL_membarrier __asm__ __volatile__ ("dmb")
/*      #define ATL_USEPCA 1 */
   #endif
#elif defined(ATL_ARCH_IA64Itan) || defined(ATL_ARCH_IA64Itan2)
   #ifdef __GNUC__
      #define ATL_membarrier __asm__ __volatile__ ("mf")
/*      #define ATL_USEPCA 1 */
   #endif
/*
 * All known x86 machines are strongly-ordered by default (can override
 * on PHI using special instructions).
 */
#elif defined(ATL_GAS_x8664) || defined (ATL_GAS_x8632)
   #define ATL_membarrier
   #define ATL_USEPCA 1
#else
   #define ATL_membarrier
#endif

#endif
