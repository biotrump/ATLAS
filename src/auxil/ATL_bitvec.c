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
#include "atlas_bitvec.h"
/*
 **************************************************************************
 * These guys used a lot in parallel programming.  When swimming through  *
 * a sea of race conditions, it can be very easy to doubt whether these   *
 * optimized routines are working correctly.  So, if you define BV_DEBUG  *
 * most of the optimized routines will check themselves against much      *
 * simpler implementations.  Debug option is nice if you want to further  *
 * optimize or add routines, so I keep it around.                         *
 * IMPORTANT: ATLAS does a lot of bit checking without mutex locks, so    *
 *    a difference in answers while running in parallel does not mean a BV*
 *    function is broken!  So, do this testing with the bitvector locked  *
 *    or in serial!                                                       *
 * #define BV_DEBUG 1                                                     *
 **************************************************************************
 */
/*
 * A bitvector consists of an integer array.  The first element of the
 * array indicates the number of full BVs required, and the 2nd elt stores
 * the remainder.
 * Assume 32 bits per int for portability
 */
int *ATL_NewBV
(
   int nmax    /* max # of bits to store */
)
{
   int *bv;
   int nelt, nrem;
   ATL_assert(sizeof(int) >= 4);
   nelt = (nmax+31) >> 5;
   nrem = nmax - ((nmax>>5)<<5);
   bv = calloc(nelt+2, sizeof(int));
   ATL_assert(bv);
   *bv = (nrem) ? nelt-1 : nelt;
   bv[1] = nrem;
/*
 * Set all bits above nrem so that last entry can be handled like all
 * others in ATL_FindFirstUnsetBitBV
 */
   if (nrem)
      bv[nelt+1] = ~((1<<nrem)-1);
   return(bv);
}


void ATL_SetBitBV
(
   int *bv,   /* bitvector to use */
   int bpos   /* which bit to set */
)
{
   int iv, i, k;
   iv = bpos >> 5;
   i = bpos - (iv<<5);
   bv[iv+2] |= 1<<i;
}

void ATL_UnsetBitBV
(
   int *bv,   /* bitvector to use */
   int bpos   /* which bit to set */
)
{
   int iv, i, k;
   iv = bpos >> 5;
   i = bpos - (iv<<5);
   bv[iv+2] &= ~(1<<i);
}

/*
 * RETURNS: bit at position bpos
 */
int ATL_IsBitSetBV
(
   int *bv,   /* bitvector to use */
   int bpos   /* which bit to set */
)
{
   int iv = bpos>>5, i = bpos-(iv<<5);
   return(bv[iv+2] & (1<<i));
}

/*
 * Returns 1 if all bits in range [b0,bn] are set, 0 otherwise
 */
#ifdef BV_DEBUG
static int ATL_IsBitRangeSetBV00
#else
int ATL_IsBitRangeSetBV
#endif
(
   int *bv,   /* bitvector to use */
   int b0,    /* first bit that must be set */
   int bN     /* last bit that must be set */
)
{
   int iv = b0>>5, i = b0-(iv<<5), nchk=0, nfvec, nr;
/*
 * False if range include bits that don't exist
 */
   if (b0 < 0 || bN < b0 || bN >= (bv[0]<<5)+bv[1])
      return(0);
    nr = bN-b0+1;   /* # of bits remaining to check */
/*
 * If we index past beginning of first bitvec, must ignore low bits.
 * No problem if partial on high end, since we store 1s in partial high end
 */
   if (i)
   {
      int v = bv[2+iv];
      int mask = (0xFFFFFFFF << i);  /* i bits at bottom are unset */
      if (nr <= 32-i) /* all bits needing to be checked are in this BV! */
      {
         int n = 32-i-nr;  /* # of bits to ignore at top */
         if (n)
            mask &= 0x7FFFFFFFF >> (n-1);
         return((v|mask) == v);
      }
/*
 *    If we get here, we are checking all upper bits in this vector
 */
      if ((v | mask) != v)
         return(0);
      nr -= 32-i;    /* we checked 32-i bits, so dec num remaining */
      if (nr < 1)
         return(1);
      iv++;          /* done checking this BV, go to next */
   }
/*
 * Remaining vectors can be compared with 0xFFFFFFFF, since we store 1s above
 * end of last partial vector
 */
   nfvec = (nr>>5);         /* # of full vectors in remaining range */
   nr -= (nfvec<<5);        /* # of bits to check in last partial vector */
/*
 * All full vectors checked wt 0xFFFFFFFF
 */
   for(i=0; i < nfvec; i++)
   {
      int v = bv[i+iv];
      if ((v|0xFFFFFFFF) != v)
         return(0);
   }
/*
 * Now, just check that last nr bits are set in any remaining vector
 */
   if (nr)
   {
      int v = bv[nfvec+iv];
      unsigned mask = ~((1<<nr)-1);  /* mask has lower nr bits set */
      if ((v|mask) != v)
         return(0);
   }
   return(1);
}
#ifdef BV_DEBUG
int ATL_IsBitRangeSetBV
(
   int *bv,   /* bitvector to use */
   int b0,    /* first bit that must be set */
   int bN     /* last bit that must be set */
)
{
   int i, iret=1;
   for (i=b0; i <= bN; i++)
   {
      if (!ATL_IsBitSetBV(bv,i))
      {
         iret = 0;
         break;
      }
   }
   if (iret)
   {
      ATL_assert(ATL_IsBitRangeSetBV00(bv, b0, bN));
   }
   else
   {
      ATL_assert(!ATL_IsBitRangeSetBV00(bv, b0, bN));
   }
}
#endif

/*
 * Returns smallest bit position that has not yet been set, -1 if all set
 * Unused bits in last bv are set to make this search fast.
 */
#ifdef BV_DEBUG
static int ATL_FindFirstUnsetBitBV00
#else
int ATL_FindFirstUnsetBitBV
#endif
(
   int *bv,   /* bitvector to use */
   int bs     /* bit to start at (0 means beginning) */
)
{
   const int nbv = bv[1] ? bv[0]+1:bv[0];
   int iv = bs>>5, is=bs-(iv<<5);
   int i;
/*
 * No match if start is beyond end of bits
 */
   if (bs >= (bv[0]<<5) + bv[1])
      return(-1);
/*
 * If first vector to check is partial
 */
   if (is)
   {
      const unsigned int v = bv[iv+2], base=iv<<5;
      if ((v|(0xFFFFFFFF<<is)) != v)
      {
         const unsigned int nb=(bv[1] && iv == nbv-1) ? bv[1]:32;
         for (i=is; i < nb; i++)
            if ((v|(1<<i)) != v)
               return(base+i);
      }
      iv++;
   }
   for (; iv < nbv; iv++)
   {
      const unsigned int v = bv[iv+2], base=iv<<5;
      int kbeg;
      if (v == 0xFFFFFFFF)   /* if all bits are set */
          continue;          /* skip to next bitvec */
/*
 *    Find start of 4-bit area where first unset bit occurs
 */
      if (v != (v|0xFF))   /* unset bit in first 8 bits */
         kbeg = (v != (v|0xF)) ? 0:4;
      else if (v != (v|0xFF00))
         kbeg = (v != (v|0xF00)) ? 8:12;
      else if (v != (v|0xFF0000))
         kbeg = (v != (v|0xF0000)) ? 16:20;
      else
         kbeg = (v != (v|0xF000000)) ? 24:28;
/*
 *    Return 1st unset bit in 4-bit area starting at kbeg: at least 1 unset!
 */
      if (v != (v|(1<<kbeg)))
         return(base+kbeg);
      else if (v != (v|(1<<(kbeg+1))))
         return(base+kbeg+1);
      else if (v != (v|(1<<(kbeg+2))))
         return(base+kbeg+2);
      else
         return(base+kbeg+3);
   }
   return(-1); /* all bits are set */
}
#ifdef BV_DEBUG
int ATL_FindFirstUnsetBitBV
(
   int *bv,   /* bitvector to use */
   int bs     /* bit to start at (0 means beginning) */
)
{
   int iret=-1, i;
   const int n = (bv[0]<<5) + bv[1];
   for (i=bs; i < n; i++)
   {
      if (!ATL_IsBitSetBV(bv, i))
      {
         iret = i;
         break;
      }
   }
   i = ATL_FindFirstUnsetBitBV00(bv, bs);
   if (i != iret)
   {
      fprintf(stderr, "BV: good=%d, bad=%d, bs=%d, n=%d\n", iret, i, bs, n);
      ATL_assert(0);
   }
}
#endif
/*
 * Find the first set bit, after skipping first bs bits
 */
#ifdef BV_DEBUG
static int ATL_FindFirstSetBitBV00
#else
int ATL_FindFirstSetBitBV
#endif
(
   int *bv,   /* bitvector to use */
   int bs     /* bit to start at (0 means beginning) */
)
{
   int iv = bs>>5, is=bs-(iv<<5);
   const int nbv = bv[0];
/*
 * No match if start is beyond end of bits
 */
   if (bs >= (bv[0]<<5) + bv[1])
      return(-1);
/*
 * If first vector to check is partial
 */
   if (is)
   {
      const unsigned int v = bv[iv+2], nchk=(iv==nbv) ? bv[1] : 32;
      if (v & ((0xFFFFFFFF)<<is))
      {
         int i;
         for (i=is; i < nchk; i++)
            if (v & (1<<i))
               return(((iv)<<5)+i);
      }
      iv++;  /* no set bits in initial bv checked */
   }
   for (; iv < nbv; iv++)
   {
      const unsigned int v = bv[iv+2], base=iv<<5;
      int kbeg;
      if (!v)                /* if no bits set */
         continue;           /* skip to next bitvec */
/*
 *    Establish start of 4-bit area where first set bit occurs
 */
      if (v&0xFF)             /* set bit in first 8 bits! */
          kbeg = (v&0xF) ? 0:4;
      else if (v&0xFF00)
          kbeg = (v&0xF00) ? 8:12;
      else if (v&0xFF0000)
          kbeg = (v&0xF0000) ? 16:20;
      else
          kbeg = (v&0xF000000) ? 24:28;
/*
 *    Ret 1st set in 4-bit area starting at kbeg where at least one bit is set
 */
      if (v & (1<<kbeg))
         return(base+kbeg);
      else if (v & 1<<(kbeg+1))
         return(base+kbeg+1);
      else if (v & 1<<(kbeg+2))
         return(base+kbeg+2);
      else
         return(base+kbeg+3);

   }
/*
 * If we have a partial bitvector at end, check it after masking off
 * unneeded bits
 */
   if (bv[1])
   {
      const unsigned int v = bv[iv+2]&(0x7FFFFFFF>>(31-bv[1])), base=iv<<5;
      int kbeg;
      if (!v)
         return(-1);
/*
 *    Establish start of 4-bit area where first set bit occurs
 */
      if (v&0xFF)             /* set bit in first 8 bits! */
          kbeg = (v&0xF) ? 0:4;
      else if (v&0xFF00)
          kbeg = (v&0xF00) ? 8:12;
      else if (v&0xFF0000)
          kbeg = (v&0xF0000) ? 16:20;
      else
          kbeg = (v&0xF000000) ? 24:28;
/*
 *    Ret 1st set in 4-bit area starting at kbeg where at least one bit is set
 */
      if (v & (1<<kbeg))
         return(base+kbeg);
      else if (v & 1<<(kbeg+1))
         return(base+kbeg+1);
      else if (v & 1<<(kbeg+2))
         return(base+kbeg+2);
      else
         return(base+kbeg+3);
   }
   return(-1);  /* no bits after bs are set */
}
#ifdef BV_DEBUG
int ATL_FindFirstSetBitBV
(
   int *bv,   /* bitvector to use */
   int bs     /* bit to start at (0 means beginning) */
)
{
   int iret=-1, i;
   const int n = (bv[0]<<5) + bv[1];
   for (i=0; i < n; i++)
   {
      if (ATL_IsBitSetBV(bv, i))
      {
         iret = i;
         break;
      }
   }
   ATL_assert(iret = ATL_FindFirstSetBitBV00(bv, bs));
}
#endif

void ATL_UnsetAllBitsBV(int *bv)
{
   const int nbv = bv[0], nr = bv[1];
   unsigned int *bp = bv+2;
   int i;
   for (i=0; i < nbv; i++)
      bp[i] = 0;
   if (nr)
      bp[nbv] = ~((1<<nr)-1);
}

void ATL_SetAllBitsBV(int *bv)
{
   const int nbv = bv[1] ? bv[0]+1 : bv[0];
   int *bp = bv+2;
   int i;
   for (i=0; i < nbv; i++)
      bp[nbv] = -1;
}

void ATL_ReverseAllBitsBV(int *bv)
{
   const int nbv = bv[0], nr = bv[1];
   unsigned int *bp = bv+2;
   int i;
   for (i=0; i < nbv; i++)
      bp[i] = ~bp[i];
   if (nr)
      bp[nbv] = (~bp[nbv]) | ((0xFFFFFFFF)<<nr);
}
