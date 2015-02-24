#include "atlas_tvec.h"

void PrintUsage(char *name, char *arg, int i)
{
   fprintf(stderr,
"This routine takes tvecs in two tvec files and produces their ratio\n");
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' ON %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s <flags> ; flags include:\n", name);
   fprintf(stderr, "   -i <file> : (stdin) file with two concatonated tvecs files to read\n");
   fprintf(stderr, "   -o <file>  : (stdout) output file for all tvecs\n");
   fprintf(stderr, "   -R # <nam1> ... <nam#>: vectors coming from both files "
           "to use for ratio\n");
   fprintf(stderr,
      "   -c # <nam1> ... <nam#>: vectors where we take first instance only\n");
   exit (i ? i : -1);
}

char **GetFlags         /* RETURNS: array of single/repeated names */
(
   int nargs,
   char **args,
   int *none,           /* # of vecs where we take only 1st definition */
   int *nmul,           /* # of vecs where we take all vectors of that name */
   FILE **fpin,         /* input stream */
   FILE **fpout         /* output stream */
)
{
   char **vc=NULL, **vr=NULL, **vv, *sp;
   int i, j, n, nc=0, nr=0;
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
      case 'c':    /* -R # <nam1> ... <nam#> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -c ", i-1);
         nr = atoi(args[i]);
         assert(nr > 0 && nr < 2048);
         vr = malloc(sizeof(char*)*nr);
         assert(vr);
         for (j=0; j < nr; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -c ", i-1);
            vr[j] = args[i];
         }
         break;
      case 'R':    /* -R # <nam1> ... <nam#> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -S ", i-1);
         nc = atoi(args[i]);
         assert(nc > 0 && nc < 2048);
         vc = malloc(sizeof(char*)*nc);
         assert(vc);
         for (j=0; j < nc; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -S ", i-1);
            vc[j] = args[i];
         }
         break;
      default :
         PrintUsage(args[0], args[i], i);
      }                                         /* end switch over flags */
   }                                            /* end for over flags */
   if (!nr && !nc)
      exit(0);   /* no output at all! */

   n = nr + nc;
   vv = malloc(n*sizeof(char*));
   assert(vv);
   for (i=0; i < nr; i++)
      vv[i] = vr[i];
   if (vr)
      free(vr);
   if (vc)
   {
      for (; i < n; i++)
         vv[i] = vc[i-nr];
      free(vc);
   }
   *none = nr;
   *nmul = nc;
   return(vv);
}

int main(int nargs, char **args)
{
   FILE *fpin, *fpout;
   char **vnams1, **vnamsr, *cmnt, *cm;
   int Nf, N1, Nr, N, i, j, RNGINC=0;
   ATL_tvec_t *ib1, *ib2, *op, *ob=NULL;

   vnams1 = GetFlags(nargs, args, &N1, &Nr, &fpin, &fpout);
   vnamsr = vnams1 + N1;
/*
 * Read in both file's tvecs as input
 */
   ib1 = ATL_ReadTvecFile(fpin, &cmnt, &N);
   ib2 = ATL_ReadTvecFile(fpin, &cm, &N);
   free(cm);
   if (fpin != stdin)
      fclose(fpin);
/*
 * Grab one copy of all -c tvecs from input, and place on output queue
 */
   for (i=0; i < N1; i++)
   {
      if (!(op=ATL_PullNamedVecsFromList(1, vnams1+i, &ib1)))
      {
          assert((op=ATL_PullNamedVecsFromList(1, vnams1+i, &ib2)));
      }
      op->next = ob;
      ob = op;
   }
/*
 * Grab one copy of all -R tvecs from each file, and create a ratio tvec
 * that is put on queue
 */
   for (i=0; i < Nr; i++)
   {
      ATL_tvec_t *op2;
      int N;
      double *d1, *d2;

      op = ATL_PullNamedVecsFromList(1, vnamsr+i, &ib1);
      op2 = ATL_PullNamedVecsFromList(1, vnamsr+i, &ib2);
      assert(op);
      assert(op2);
      assert(op->pre == 'd' && op2->pre == 'd');
      assert(op->nrep = op2->nrep);
      assert(op->N == op->N);
      N = op->N;
      d1 = op->vp;
      d2 = op2->vp;
      for (j=0; j < N; j++)
      {
         d1[j] /= d2[j];
      }
      ATL_KillThisTvec(op2);
      op->next = ob;
      ob = op;
   }
   ATL_KillAllTvecs(ib1);
   ATL_KillAllTvecs(ib2);
/*
 * Get list in order from input: -c tvecs -S tvecs
 */
   ob = ATL_ReverseTvecList(ob);
/*
 * Write them out, and we are done
 */
   ATL_WriteTvecFile(fpout, cmnt, ATL_CountTvecsInList(ob), ob);
   ATL_KillAllTvecs(ob);
   free(cmnt);
   free(vnams1);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
   return(0);
}
