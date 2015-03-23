/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2014 Md Rakib Hasan
 *
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
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
   fprintf(stderr, "This routine takes tvecs in one or more files and outputs"
                   " the \nselected ones, while optionally renaming them.\n");
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' ON %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s <flags> ; flags include:\n", name);
   fprintf(stderr, "   -i <file> : (stdin) file with vecs to read.\n");
   fprintf(stderr, "   -o <file> : (stdout) output file for all tvecs.\n");
   fprintf(stderr, "   -# <#> : # of tvec files concatonated in the input; "
                   "MUST be specified\n");
   fprintf(stderr, "      before -S or -R.\n");
   fprintf(stderr, "   -S # <nam1> ... <nam#> : vectors to select from 1st "
                   "file.\n      Repeat for later files.\n");
   fprintf(stderr, "   -R # <oldnam1> <newnam1> ... <oldnam#> <newnam#> : "
                   "name pairs from 1st file.\n      Used with -S. Repeat flag "
                   "for later files.\n");
   fprintf(stderr, "   -s # <nam1> ... <nam#> : vectors to select from all "
                   "files.\n      Can be used once. Cannot be used with "
                   "-S or -R.\n");
   fprintf(stderr, "   -r # <oldnam1> <newnam1> ... <oldnam#> <newnam#> : "
                   " rename list\n      for all files. "
                   "Can be used once and only with -s.\n");
   fprintf(stderr, "\n");
   exit (i ? i : -1);
}

void GetFlags
(
   int nargs,
   char **args,
   int *Nfiles,         /* # of tvec files in input stream */
   char *Mode,          /* S or s for -S/R or -s/r */
   int **Nselect,       /* array holding # of tvecs to select for each file */
   char ****sel_names,  /* array holding names of tvecs for each file */
   int **Nrename,       /* array holding # of tvecs to rename for each file */
   char ****ren_names,  /* array holding names of tvecs for each file */
   FILE **fpin,         /* input stream */
   FILE **fpout         /* output stream */
)
{
   int i, j, k, n;
   int si, ri;

   si = ri = 0;
   *fpin = stdin;
   *fpout = stdout;
   *Nfiles = 1;
   *Mode = 0; /* initially no mode */
   *Nselect = malloc(*Nfiles * sizeof(int));
   *sel_names = malloc(*Nfiles * sizeof(char**));
   *Nrename = malloc(*Nfiles * sizeof(int));
   *ren_names = malloc(*Nfiles * sizeof(char**));
   for (k=0; k<*Nfiles; k++)
   {
      (*Nselect)[k] = 0;
      (*Nrename)[k] = 0;
   }
   for (i=1; i<nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], "no '-' preceeding flag!", i);

      switch(args[i][1])
      {
         case '#':   /* -# <# files> */
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -# ", i-1);
            *Nfiles = atoi(args[i]);
            free(*Nselect);
            free(*sel_names);
            free(*Nrename);
            free(*ren_names);
            *Nselect = malloc(*Nfiles * sizeof(int));
            *sel_names = malloc(*Nfiles * sizeof(char**));
            *Nrename = malloc(*Nfiles * sizeof(int));
            *ren_names = malloc(*Nfiles * sizeof(char**));
            for (k=0; k<*Nfiles; k++)
            {
               (*Nselect)[k] = 0;
               (*Nrename)[k] = 0;
            }
            break;
         case 'i':   /* -i <file> */
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -i ", i-1);
            *fpin = fopen(args[i], "r");
            assert(*fpin);
            break;
         case 'o':   /* -o <file> */
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -o ", i-1);
            *fpout = fopen(args[i], "w");
            assert(*fpout);
            break;
         case 'S':   /* -S # <nam1> ... <nam#> */
         case 's':   /* -s # <nam1> ... <nam#> */
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -S/s ", i-1);
            assert(*Mode == 0 || *Mode == args[i-1][1]);
            *Mode = args[i-1][1];
            if (*Mode == 's') assert(si == 0);
            (*Nselect)[si] = n = atoi(args[i]);
            if (n > 0)
            {
               (*sel_names)[si] = malloc(n * sizeof(char*));
               assert((*sel_names)[si]);
               for (j=0; j<n; j++)
               {
                  if (++i >= nargs)
                     PrintUsage(args[0], "out of flags in -S/s ", i-1);
                  (*sel_names)[si][j] = args[i];
               }
            }
            si++;
            break;
         case 'R':   /* -R # <oldnam1> <newnam1> ... <oldnam#> <newnam#> */
         case 'r':   /* -r # <oldnam1> <newnam1> ... <oldnam#> <newnam#> */
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -R/r ", i-1);
            assert(*Mode == 0 || *Mode == (args[i-1][1]+1));
            *Mode = args[i-1][1]+1;  /* mode s/S */
            if (*Mode == 'r') assert(ri == 0);
            (*Nrename)[ri] = n = atoi(args[i]);
            if (n > 0)
            {
               (*ren_names)[ri] = malloc(n*2 * sizeof(char*));
               assert((*ren_names)[ri]);
               for (j=0; j<n*2; j++) /* pair of names */
               {
                  if (++i >= nargs)
                     PrintUsage(args[0], "out of flags in -R/r ", i-1);
                  (*ren_names)[ri][j] = args[i];
               }
            }
            ri++;
            break;
         default:
            PrintUsage(args[0], args[i], i);
      }
   }
}

int main(int nargs, char **args)
{
   FILE *fpin, *fpout;
   int Nf, *Ns=NULL, *Nr=NULL;
   int i, j, n, ci;
   char mode = 0;
   char ***snames=NULL, ***rnames=NULL, *cmnt0=NULL, *cmnt=NULL;
   ATL_tvec_t *tp, *np, *nb, *ft;

   GetFlags(nargs, args, &Nf, &mode, &Ns, &snames, &Nr, &rnames, &fpin, &fpout);

   if (Ns[0] > 0)
   {
      tp = ATL_ReadTvecFile(fpin, &cmnt0, &n);
      /* select the tvecs */
      nb = ATL_PullNamedVecsFromListWithDups(Ns[0], snames[0], &tp);
      assert(nb);
      if (tp) ATL_KillAllTvecs(tp);
      /* rename the tvecs */
      for (j=0; j<Nr[0]; j++)
      {
         tp = ATL_FindTvecByName(nb, rnames[0][j*2]);
         if (tp) strcpy(tp->name, rnames[0][j*2+1]);
      }
   }
   for (i=1; i<Nf; i++)
   {
      ci = ((mode == 'S') ? i : 0);
      if (Ns[ci] > 0)
      {
         tp = ATL_ReadTvecFile(fpin, &cmnt, &n);
         np = ATL_PullNamedVecsFromListWithDups(Ns[ci], snames[ci], &tp);
         if (!np) break;
         if (tp) ATL_KillAllTvecs(tp);
         if (cmnt)
         {
            free(cmnt);
            cmnt = NULL;
         }

         /* rename the tvecs */
         for (j=0; j<Nr[ci]; j++)
         {
            tp = ATL_FindTvecByName(np, rnames[ci][j*2]);
            if (tp) strcpy(tp->name, rnames[ci][j*2+1]);
         }

         /* append the current list to the base list */
         ATL_FindLastTvecInList(nb)->next = np;
      }
   }
   if (mode == 's') /* -s/r mode, remove duplicates */
   {
      for (j=0; j<Nr[0]; j++)
      {
         for (i=0; i<Ns[0]; i++)
         {
            if (!strcmp(rnames[0][j*2], snames[0][i]))
            {
               strcpy(snames[0][i], rnames[0][j*2+1]);
            }
         }
      }
      tp = ATL_PullNamedVecsFromList(Ns[0], snames[0], &nb);
      ATL_KillAllTvecs(nb);
      nb = tp;
   }

   /*
    * Write them out, and we are done
    */
   ATL_WriteTvecFile(fpout, cmnt0, ATL_CountTvecsInList(nb), nb);
   ATL_KillAllTvecs(nb);

   /* free the args */
   for (i=0; i<Nf; i++)
   {
      if (Ns && Ns[i] > 0) free(snames[i]);
      if (Nr && Nr[i] > 0) free(rnames[i]);
   }
   if (cmnt0) free(cmnt0);
   if (snames) free(snames);
   if (rnames) free(rnames);
   if (Ns) free(Ns);
   if (Nr) free(Nr);
   if (fpin != stdin) fclose(fpin);
   if (fpout != NULL && fpout != stdout && fpout != stderr)
      fclose(fpout);
   return(0);
}
