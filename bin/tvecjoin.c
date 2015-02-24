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
#include "atlas_tvec.h"

void PrintUsage(char *name, char *arg, int i)
{
   fprintf(stderr, "Joins multiple vectors so they form a contiguous range.\n");
   fprintf(stderr, "Range vectors should not overlap.\n\n");
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' ON %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s <flags> ; flags include:\n", name);
   fprintf(stderr, "   -i <file> : (stdin) file wt vecs to combine\n");
   fprintf(stderr, "   -o <file>  : (stdout) file for joined vecs\n");
   fprintf(stderr, "   -J # <rngv1> <joinv1> ... <rngv#> <joinV#>\n");
   fprintf(stderr,
      "      Vectors are joined into the first -J pair (rangv1 joinv1)\n");
   fprintf(stderr,
           "      The -J flag may be repeated as many times as necessary.\n");
   exit (i ? i : -1);
}

typedef struct ATL_tvjoin ATL_tvjoin_t;
struct ATL_tvjoin
{
   int N;
   char **rangvs, **joinvs;
   ATL_tvjoin_t *next;
};

ATL_tvjoin_t *ATL_tvj_new(int N)
{
   ATL_tvjoin_t *tp;
   tp = malloc(sizeof(ATL_tvjoin_t));
   assert(tp);
   tp->N = N;
   tp->rangvs = calloc(N+N,sizeof(char *));
   assert(tp->rangvs);
   tp->joinvs = tp->rangvs + N;
   return(tp);
}

ATL_tvjoin_t *ATL_tvj_kill(ATL_tvjoin_t *tp)
{
   ATL_tvjoin_t *tn=NULL;
   if (tp)
   {
      tn = tp->next;
      free(tp->rangvs);  /* these pt to commandline args, so don't free them */
      free(tp);
   }
   return(tn);
}

ATL_tvjoin_t *GetFlags   /* RETURNS: queue of joinings to perform */
(
   int nargs,
   char **args,
   FILE **fpin,         /* input stream */
   FILE **fpout         /* output stream */
)
{
   char **na=NULL, *sp;
   int i, j, n;
   FILE *fp;
   ATL_tvjoin_t *jb=NULL, *jp, *jr;

   *fpin = stdin;
   *fpout = stdout;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], "no '-' preceeding flag!", i);
      switch(args[i][1])
      {
      case 'i':    /* -i[1,2] <file> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -i ", i-1);
         fp = fopen(args[i], "r");
         assert(fp);
         *fpin = fp;
         break;
      case 'o':    /* -o <file> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -i ", i-1);
         fp = fopen(args[i], "w");
         assert(fp);
         *fpout = fp;
         break;
      case 'J':  /* -J # <rngv1> <joinv1> ... <rngv#> <joinV#> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -J ", i-1);
         n = atoi(args[i]);
         assert(n > 1 && n < 2048);
         jp = ATL_tvj_new(n);
         for (j=0; j < n; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -J ", i-1);
            jp->rangvs[j] = args[i];
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -J ", i-1);
            jp->joinvs[j] = args[i];
         }
         jp->next = jb;
         jb = jp;
         break;
      default :
         PrintUsage(args[0], args[i], i);
      }                                         /* end switch over flags */
   }                                            /* end for over flags */
   if (!jb)
   {
      fprintf(stderr, "You must join something!\n");
      exit(-1);
   }
/*
 * Reverse the join order so its in commandline order!
 */
   jr = NULL;
   while (jb)
   {
      jp = jb->next;
      jb->next = jr;
      jr = jb;
      jb = jp;
   }
   return(jr);
}

int main(int nargs, char **args)
{
   char *cmnt;
   FILE *fpin, *fpout;
   ATL_tvjoin_t *jb, *jp;
   ATL_tvec_t *tbu, *tjb=NULL;
   int i;

   jb = GetFlags(nargs, args, &fpin, &fpout);
/*
 * Read in all unjoined vectors
 */
   tbu = ATL_ReadTvecFile(fpin, &cmnt, &i);
/*
 * Loop over all joined vectors that need to be created
 */
   for (jp=jb; jp; jp = jp->next)
   {
      ATL_tvec_t *tb, *tp, *rb, *td, *rd;
      const int N = jp->N;
      int j;
/*
 *    Yank join & dup range pointers from unjoined pool
 */
      rb = ATL_DupNamedVecsFromList(N, jp->rangvs, tbu, 1);
      assert(ATL_CountTvecsInList(rb) == N);
      tb = ATL_PullNamedVecsFromList(N, jp->joinvs, &tbu);
      assert(ATL_CountTvecsInList(tb) == N);
/*
 *    Yank destination registers from list
 */
      rd = ATL_PullNamedVecsFromList(1, jp->rangvs, &rb);
      td = ATL_PullNamedVecsFromList(1, jp->joinvs, &tb);
      assert(rd->pre == 'i');
/*
 *    Loop over source registers, and join them into destination registers
 */
      for (j=1; j < N; j++)
      {
         ATL_tvec_t *ts, *rs;  /* source vectors */
/*
 *       Find source range & tvec, then join them into destination
 */
         rs = ATL_FindTvecByName(rb, jp->rangvs[j]);
         ts = ATL_FindTvecByName(tb, jp->joinvs[j]);
         tp = ATL_CombineTheseVecsUsingInts(rd, rs, td, ts);
         ATL_KillThisTvec(td);
         td = tp;
         tp = ATL_CombineTheseVecsUsingInts(rd, rs, rd, rs);
         ATL_KillThisTvec(rd);
         rd = tp;
      }
/*
 *    Kill the original tvecs that have been joined
 */
      ATL_KillAllTvecs(tb);
      ATL_KillAllTvecs(rb);
/*
 *    Add joined tvecs to queue.  Range vectors may be repeated, so don't
 *    add if it is already there.
 */
      td->next = tjb;
      tjb = td;
      tp = ATL_FindTvecByName(tjb, rd->name);
      if (!tp)
      {
         rd->next = tjb;
         tjb = rd;
      }
      else /* check that usage is correct */
      {
         int *ip = rd->vp, *ip2=tp->vp;
         const int N = rd->N;
         assert(tp->N == N && tp->nrep == rd->nrep && tp->pre == rd->pre);
         for (i=0; i < N; i++)
            assert(ip[i] == ip2[i]);
      }
   }

   ATL_KillAllTvecs(tbu);
   ATL_WriteTvecFile(fpout, cmnt, ATL_CountTvecsInList(tjb), tjb);
   ATL_KillAllTvecs(tjb);
   free(cmnt);
   while (jb)
      jb = ATL_tvj_kill(jb);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
   return(0);
}
