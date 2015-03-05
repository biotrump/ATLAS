#include "atlas_tvec.h"

void PrintUsage(char *name, char *arg, int i)
{
   fprintf(stderr,
"This routine takes tvecs in multiple files and puts then into one\n");
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' ON %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s <flags> ; flags include:\n", name);
   fprintf(stderr, "   -i <file> : (stdin) file with vecs to read\n");
   fprintf(stderr, "   -o <file>  : (stdout) output file for all tvecs\n");
   fprintf(stderr,
      "   -C # <nam1> ... <nam#>: vectors coming from all files\n");
   fprintf(stderr,
      "   -c # <nam1> ... <nam#>: vectors where we take first instance only\n");
   fprintf(stderr,
      "   -f # <nam1> ... <nam#>: vectors where we take only\n"
      "                           the first instance from each file\n");
   fprintf(stderr, "   -s 0/1 : don't/do suffix tvecs wt file number\n");
   fprintf(stderr,
           "   -S # suf1 ... sufN: suffix first N files' tvecs wt "
           "provided strings.\n"
           "                       Overrides -s for first N files.\n");
   exit (i ? i : -1);
}

ATL_tvstrq_t *GetFlags         /* RETURNS: queue of file suffixes to use */
(
   int nargs,
   char **args,
   int *SUFFCNT,
   ATL_tvstrq_t **namr, /* names to take all repeats (-C) */
   ATL_tvstrq_t **nam1, /* names to take first instance only (-c) */
   ATL_tvstrq_t **namf, /* names to take first instance frm each file (-f) */
   FILE **fpin,         /* input stream */
   FILE **fpout         /* output stream */
)
{
   char **vc=NULL, **vr=NULL, *vf, **vv, *sp;
   int i, j, n, nc=0, nr=0;
   FILE *fp;
   ATL_tvstrq_t *qp, *suffb=NULL;

   *SUFFCNT=0;
   *namr = *nam1 = *namf = NULL;
   *fpin = stdin;
   *fpout = stdout;
   for (i=1; i < nargs; i++)
   {
      int SKIP=0;
      ATL_tvstrq_t **qb=NULL;
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
      case 'C':    /* -c # <nam1> ... <nam#> */
         qb = namr;
      case 'f':    /* -f # <nam1> ... <nam#> */
         if (!qb)
            qb = namf;
      case 'c':    /* -c # <nam1> ... <nam#> */
         if (!qb)
            qb = nam1;
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -c ", i-1);
         if (isdigit(args[i][0]))
            nr = atoi(args[i]);
         else
            SKIP = nr = 1;
         assert(nr > 0 && nr < 2048);
         for (j=0; j < nr; j++)
         {
            if (!SKIP)
            {
               if (++i >= nargs)
                  PrintUsage(args[0], "out of flags in -c ", i-1);
               SKIP=0;
            }
            qp = ATL_GetStrNode(args[i]);
            qp->next = *qb;
            *qb = qp;
         }
         break;
      case 'S':    /* -S # <suff1> ... <suff#> */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -c ", i-1);
         if (isdigit(args[i][0]))
            nr = atoi(args[i]);
         else
            SKIP = nr = 1;
         assert(nr > 0 && nr < 2048);
         for (j=0; j < nr; j++)
         {
            if (!SKIP)
            {
               if (++i >= nargs)
                  PrintUsage(args[0], "out of flags in -c ", i-1);
               SKIP=0;
            }
            qp = ATL_GetStrNode(args[i]);
            qp->next = suffb;
            suffb = qp;
         }
         break;
      case 's': /* -s [0,1] */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -s ", i-1);
         *SUFFCNT = (args[i][0] != '0');
         break;
      default :
         PrintUsage(args[0], args[i], i);
      }                                         /* end switch over flags */
   }                                            /* end for over flags */
   if (suffb)
      suffb = ATL_ReverseStrQ(suffb);
   if (*namf)
      *namf = ATL_ReverseStrQ(*namf);
   if (*nam1)
      *nam1 = ATL_ReverseStrQ(*nam1);
   if (*namr)
      *namr = ATL_ReverseStrQ(*namr);

   return(suffb);
}

int main(int nargs, char **args)
{
   FILE *fpin, *fpout;
   char *cmnt;
   ATL_tvec_t *tb;
   ATL_tvstrq_t *sfb, *nrb, *n1b, *nfb;
   int SUFFCNT;
   ATL_tvec_t *ob=NULL;

   sfb = GetFlags(nargs, args, &SUFFCNT, &nrb, &n1b, &nfb, &fpin, &fpout);
/*
 * If I'm not selecting tvecs, then I'll be outputing all tvecs
 */
   if (!nrb && !n1b && !nfb)
   {
      int nF=0, i;
      while ( (tb = ATL_ReadTvecFile(fpin, &cmnt, &i)) )
      {
         nF++;
         free(cmnt);
         if (sfb)
         {
            ATL_SuffixTvecNames(tb, sfb->sp);
            sfb = ATL_KillStrNode(sfb);
         }
         else if (SUFFCNT)
         {
            char suff[32];
            sprintf(suff, "_%d", nF);
            ATL_SuffixTvecNames(tb, suff);
         }
         if (ob)
            ATL_FindLastTvecInList(ob)->next = tb;
         else
            ob = tb;
      }
   }
   else  /* I'm pulling out only selected tvecs */
   {
      char **vnams1=NULL, **vnamsR=NULL, *cmnt;
      int nR=0, nf1=0, nF=0, i;
      ATL_tvec_t *b=NULL;

      if (nrb)
         vnamsR = ATL_StrQ2Arr(nrb, &nR);
      if (nfb)
         vnams1 = ATL_StrQ2Arr(nfb, &nf1);
      while ( (tb = ATL_ReadTvecFile(fpin, &cmnt, &i)) )
      {
         ATL_tvec_t *pp, *b=NULL;  /* pull ptr */
         free(cmnt);
         nF++;
/*
 *       See if there are any 1-time vectors to take
 */
         if (n1b)
         {
/*
 *          Handle cases that match on n1b
 */
            while ( (pp = ATL_FindTvecByName(tb, n1b->sp)) )
            {
               tb = ATL_RemoveTvecFromList(tb, pp);
               pp->next = NULL;
               if (b)
                  ATL_FindLastTvecInList(b)->next = pp;
               else
                  b = pp;
               n1b = ATL_KillStrNode(n1b);
               if (!n1b)
                  break;
            }
/*
 *          Process non-base 1-time tvecs; removing them does not affect n1b
 */
            if (n1b)
            {
               ATL_tvstrq_t *qp=n1b->next;
               while (qp)
               {
                  pp = ATL_FindTvecByName(tb, qp->sp);
                  if (pp)
                  {
                     tb = ATL_RemoveTvecFromList(tb, pp);
                     pp->next = NULL;
                     if (b)
                        ATL_FindLastTvecInList(b)->next = pp;
                     else
                        b = pp;
                     qp = ATL_KillStrNode(qp);
                  }
                  else    /* didn't find tvec I'm looking for */
                     qp = qp->next;
               }
            }    /* end n1b if */
         }       /* end of if for 1-time tvec removal */
         if (nR)
         {
            pp = ATL_PullNamedVecsFromListWithDups(nR, vnamsR, &tb);
            if (b)
               ATL_FindLastTvecInList(b)->next = pp;
            else
               b = pp;
         }
         if (nf1)
         {
            pp = ATL_PullNamedVecsFromList(nf1, vnams1, &tb);
            if (b)
               ATL_FindLastTvecInList(b)->next = pp;
            else
               b = pp;
         }
/*
 *       At this point, we have removed all the required tvecs from the
 *       original base, and b has the complete tvec queue from this file
 *       in the correct order.  We now want to perform any required renaming,
 *       and then add it to the global output tvec list (ob)
 */
         if (sfb)  /* we've got user-provided suffixes to use */
         {
            ATL_SuffixTvecNames(b, sfb->sp);
            sfb = ATL_KillStrNode(sfb);
         }
         else if (SUFFCNT)  /* we need to add _file# suffix */
         {
            char suff[32];
            sprintf(suff, "_%d", nF);
            ATL_SuffixTvecNames(b, suff);
         }
         if (ob)
            ATL_FindLastTvecInList(ob)->next = b;
         else
            ob = b;
         if (tb)
            ATL_KillAllTvecs(tb);  /* get rid of this files leftover tvecs */
      } /* end while on reading in tvec files */
      if (vnams1)
      {
         for (i=0; i < nf1; i++)
            if (vnams1[i])
               free(vnams1[i]);
         free(vnams1);
      }
      if (vnamsR)
      {
         for (i=0; i < nR; i++)
            if (vnamsR[i])
               free(vnamsR[i]);
         free(vnamsR);
      }
   }    /* end if on selected or all tvecs */
   if (ob)
      ATL_WriteTvecFile(fpout, "file created by tveccat",
                        ATL_CountTvecsInList(ob), ob);
   if (fpin != stdin)
      fclose(fpin);
   ATL_KillAllTvecs(ob);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
   return(0);
}
