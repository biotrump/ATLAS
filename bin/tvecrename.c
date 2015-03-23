#include "atlas_tvec.h"

void PrintUsage(char *name, char *arg, int i)
{
   fprintf(stderr,
      "This routine renames tvecs given by -R; other vecs left alone\n");
   fprintf(stderr,
      "If -s is used, all tvecs (incl -R) are suffixed\n");
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' ON %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s <flags> ; flags include:\n", name);
   fprintf(stderr, "   -i <file> : (stdin) input file\n");
   fprintf(stderr, "   -o <file>  : (stdout) output file\n");
   fprintf(stderr, "   -s <suff> : append suffix to all names in file\n");
   fprintf(stderr,
           "   -R # <nam1> <rep1>... <nam#> <rep#>: orig & new names\n");
   exit (i ? i : -1);
}

char **GetFlags         /* RETURNS: array of orig & new names */
(
   int nargs,
   char **args,
   int *nrname,         /* # of vecs to rename */
   char **suff,
   FILE **fpin,         /* input stream */
   FILE **fpout         /* output stream */
)
{
   char **vr=NULL, **vv, *sp;
   int i, j, n, nk=0, nr=0;
   FILE *fp;

   *fpin = stdin;
   *fpout = stdout;
   *suff = NULL;
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
      case 's': /* os <suffix> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -s ", i-1);
         *suff = args[i];
         break;
      case 'R':    /* -R # <nam1> <rep1> ... <nam#> <rep#> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -R ", i-1);
         nr = atoi(args[i]);
         assert(nr > 0 && nr < 2048);
         vr = malloc(sizeof(char*)*(nr+nr));
         assert(vr);
         for (j=0; j < nr; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -R ", i-1);
            vr[j] = args[i];
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -R ", i-1);
            vr[j+nr] = args[i];
         }
         break;
      default :
         PrintUsage(args[0], args[i], i);
      }                                         /* end switch over flags */
   }                                            /* end for over flags */
   if (nr < 1 && *suff == NULL)
   {
      fprintf(stderr, "Must rename at least one vector!\n");
      exit(-1);
   }
   *nrname = nr;
   return(vr);
}

int main(int nargs, char **args)
{
   FILE *fpin, *fpout;
   char **renarr, *cmnt, *suff;
   int N, Nr, nrep, i, sL=0;
   ATL_tvec_t *tp, *tb;

   renarr = GetFlags(nargs, args, &Nr, &suff, &fpin, &fpout);
   if (suff)
      sL = strlen(suff);
/*
 * Read in all vectors in file
 */
   tb = ATL_ReadTvecFile(fpin, &cmnt, &N);
   if (fpin != stdin)
      fclose(fpin);
/*
 * Rename all targeted vectors
 */
   for (tp=tb; tp; tp = tp->next)
   {
      for (i=0; i < Nr; i++)
      {
         char *sp = renarr[i];
         if (sp && !strcmp(sp, tp->name))
         {
            int L;
            L = strlen(renarr[Nr+i])+1;
            if (tp->name)
               free(tp->name);
            tp->name = malloc(L);
            assert(tp->name);
            strcpy(tp->name, renarr[Nr+i]);
            renarr[i] = NULL;
         }
      }
      if (sL)
      {
         int L;
         char *sp = tp->name;
         L = strlen(sp);
         tp->name = malloc(L+sL+1);
         assert(tp->name);
         strcpy(tp->name, sp);
         strcpy(tp->name+L, suff);
         free(sp);
      }
   }
/*
 * Write the renamed vectors out, and we're done
 */
   ATL_WriteTvecFile(fpout, cmnt, N, tb);
   ATL_KillAllTvecs(tb);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
   if (cmnt)
      free(cmnt);
   free(renarr);
   return(0);
}
