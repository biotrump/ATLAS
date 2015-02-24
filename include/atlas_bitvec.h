/*
 *             Automatically Tuned Linear Algebra Software v3.11.31
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
#ifndef ATLAS_BITVEC_H
   #define ATLAS_BITVEC_H
   #include "atlas_misc.h"
/*
 * A bitvector consists of an integer array.  The first element of the
 * array indicates the number of full BVs required, and the 2nd elt stores the
 * number of the remainder.
 * Assumes 32 bits per int for portability
 */
#define ATL_FreeBV(bv_) free(bv_);  /* bitvec just an int array */

int *ATL_NewBV(int nmax);
void ATL_SetBitBV(int *bv, int bpos);
void ATL_UnsetBitBV(int *bv, int bpos);
int ATL_IsBitSetBV(int *bv, int bpos);
int ATL_FindFirstUnsetBitBV(int *bv, int bs);
int ATL_FindFirstSetBitBV(int *bv, int bs);
void ATL_ReverseAllBitsBV(int *bv);
void ATL_SetAllBitsBV(int *bv);
void ATL_UnsetAllBitsBV(int *bv);
int ATL_IsBitRangeSetBV(int *bv, int b0, int bN);

#endif
