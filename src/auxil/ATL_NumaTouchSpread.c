/*
 *             Automatically Tuned Linear Algebra Software v3.11.31
 *                    (C) Copyright 1999 R. Clint Whaley
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
#include "atlas_misc.h"

#include "atlas_dsysinfo.h"
#include "atlas_threads.h"
typedef struct
{
   size_t N;    /* total size in bytes */
   char *buff;  /* buffer to be spread over cores */
} ATL_TNUMA_t;

void ATL_NumaCoreTouch(ATL_LAUNCHSTRUCT_t *lp, void *vp)
{
   ATL_thread_t *tp=vp;
   ATL_TNUMA_t *np=(ATL_TNUMA_t*)lp->opstruct;
   const size_t inc = ((size_t)ATL_pgsz)*tp->P, N = np->N;
   char *cp = np->buff;
   size_t J;

   for (J=((size_t)ATL_pgsz)*tp->rank; J < N; J += inc)
      cp[J] = 0;
}

/*
 * This routine takes recently allocated buffer b, splits it into pgsz
 * chunks, and spreads the pages over all cores evenly assuming a first-touch
 * allocation strategy.  This is used to avoid having all memory owned by
 * one core, which is a disaster on a system that doesn't use a distributed
 * directory (particularly AMD HTassist).
 */
void ATL_NumaTouchSpread(size_t N, void *buff)
{
   ATL_TNUMA_t num;
   num.N = N;
   num.buff = buff;
   ATL_goparallel(ATL_NTHREADS, ATL_NumaCoreTouch, &num, NULL);
}
