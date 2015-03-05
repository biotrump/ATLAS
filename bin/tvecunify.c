#include "atlas_tvec.h"

void PrintUsage(char *name, char *arg, int i)
{
   fprintf(stderr,
      "This routine combines all tvecs with duplicate names into one\n");
   fprintf(stderr,
      "by simply stacking their entries together in the order they appear\n");
   fprintf(stderr,
      "in the input (reverse order, if -r 1 is thrown)\n");

   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' ON %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s <flags> ; flags include:\n", name);
   fprintf(stderr, "   -i <file> : (stdin) input file\n");
   fprintf(stderr, "   -o <file> : (stdout) output file\n");
   fprintf(stderr,
      "   -r 0/1    : don't / do reverse input order before stacking\n");
   exit (i ? i : -1);
}

int GetFlags         /* true if input should be reversed */
(
   int nargs,
   char **args,
   FILE **fpin,         /* input stream */
   FILE **fpout         /* output stream */
)
{
   int i, REV=0;
   FILE *fp;

   *fpin = stdin;
   *fpout = stdout;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], "no '-' preceeding flag!", i);
      switch(args[i][1])
      {
      case 'i':    /* -i <file> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -i ", i-1);
         *fpin = fopen(args[i], "r");
         assert(*fpin);
         break;
      case 'o':    /* -o <file> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -i ", i-1);
         fp = fopen(args[i], "w");
         assert(fp);
         *fpout = fp;
         break;
      case 'r': /* os <suffix> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -r ", i-1);
         REV = atoi(args[i]);
         break;
      default :
         PrintUsage(args[0], args[i], i);
      }                                         /* end switch over flags */
   }                                            /* end for over flags */
   return(REV);
}

int main(int nargs, char **args)
{
   FILE *fpin, *fpout;
   char *cmnt;
   int N, rev;
   ATL_tvec_t *tp, *tb;

   rev = GetFlags(nargs, args, &fpin, &fpout);
/*
 * Read in all vectors in file
 */
   tb = ATL_ReadTvecFile(fpin, &cmnt, &N);
   if (fpin != stdin)
      fclose(fpin);
   if (tb)
   {
      if (rev)
         tb = ATL_ReverseTvecList(tb);
      for (tp=tb; tp->next; tp = tp->next)
      {
         ATL_tvec_t *tk;
         while((tk = ATL_FindTvecByName(tp->next, tp->name)))
         {
            ATL_AdhereTvec(tp, tk);
            ATL_RemoveTvecFromList(tp, tk);
            ATL_KillThisTvec(tk);
         }
         if (!tp->next)
            break;
      }
   }
/*
 * Write the remaining vectors out, and we're done
 */
   ATL_WriteTvecFile(fpout, cmnt, ATL_CountTvecsInList(tb), tb);
   ATL_KillAllTvecs(tb);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
   if (cmnt)
      free(cmnt);
   return(0);
}
