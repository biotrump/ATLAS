/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2015 R. Clint Whaley
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
#include "atlas_bitvec.h"
void ATL_print1dBV(int col, int N, int *bv)
{
   int i, ierr;
   printf("\n\nBITVEC MAP:\n");
   ierr = ATL_FindFirstSetBitBV(bv, 0);
   if (ierr == -1)
      ierr = N;
   if (col)
   {
      for (i=0; i < N; i++)
      {
         if (i != ierr)
            printf("%d .\n", i%10);
         else
         {
            printf("%d X\n", i%10);
            if (i < N-1)
               ierr = ATL_FindFirstSetBitBV(bv, i+1);
         }
      }
   }
   else
   {
      for (i=0; i < N; i++)
         printf("%d", i%10);
      printf("\n");
      for (i=0; i < N; i++)
      {
         if (i != ierr)
            printf(".");
         else
         {
            printf("X");
            if (i < N-1)
               ierr = ATL_FindFirstSetBitBV(bv, i+1);
         }
      }
      printf("\n");
   }
}
