#ifndef ATL_Tvec_H
   #define ATL_Tvec_H

#include "atlas_genparse.h"


typedef struct ATL_Tvec ATL_tvec_t;
struct ATL_Tvec
{
   void *vp;            /* pointer to vector, type in pre */
   char *name;          /* name of vector (eg, "MFLOP" or "N") */
   ATL_tvec_t *next;
   int  N;              /* # of elements in vector (inclucing reps) */
   int nrep;            /* # of repititions in timings */
   char pre;            /* double (d), string (s), char (c), integer (i) */
};

typedef struct ATL_tvstrq ATL_tvstrq_t;
struct ATL_tvstrq
{
   char *sp;
   ATL_tvstrq_t *next;
   int len;
};

/* procedure 1 */
ATL_tvstrq_t *ATL_GetStrNode(char *name)
{
   ATL_tvstrq_t *np;
   np = malloc(sizeof(ATL_tvstrq_t));
   assert(np);
   np->sp = DupString(name);
   assert(np->sp);
   np->len = strlen(name);
   strncpy(np->sp, name, np->len+1);
   np->next = NULL;
   return(np);
}

/* procedure 2 */
ATL_tvstrq_t *ATL_KillStrNode(ATL_tvstrq_t *kp)
{
   ATL_tvstrq_t *kn=NULL;
   if (kp)
   {
      kn = kp->next;
      if (kp->sp)
         free(kp->sp);
      free(kp);
   }
   return(kn);
}

/* procedure 3 */
void ATL_KillStrQ(ATL_tvstrq_t *kb)
{
   while (kb)
      kb = ATL_KillStrNode(kb);
}

/* procedure 4 */
ATL_tvstrq_t *RemoveStrNodeFromQ(ATL_tvstrq_t *b, ATL_tvstrq_t *p)
/*
 * Safely removes p from queue b, returns possibly changed b
 */
{
   if (!b)
      return(NULL);
   if (!p)
      return(b);
   if (p == b)
      b = b->next;
   else
   {
      ATL_tvstrq_t *prev=b, *qp;
      for (qp=b->next; qp; qp = qp->next)
         prev=qp;
      prev->next = p->next;
   }
   p->next = NULL;
   return(b);
}

/* procedure 5 */
char **ATL_StrQ2Arr(ATL_tvstrq_t *qb, int *N)
/*
 * translates qb to a NULL-term array of strings, destroying qb in process
 */
{
   int n, i;
   char **sarr;
   ATL_tvstrq_t *qp;

   n = *N = ATL_CountStrNodes(qb);
   sarr = malloc((n+1)*sizeof(char*));
   assert(sarr);

   sarr[n] = NULL;
   for (i=0,qp=qb; qp; i++,qp = qp->next)
   {
      sarr[i] = qp->sp;
      qp->sp = NULL;
   }
   ATL_KillStrQ(qb);
   return(sarr);
}


/* procedure 6 */
int ATL_CountStrNodes(ATL_tvstrq_t *qb)
{
   int i;
   ATL_tvstrq_t *qp;
   for (i=0,qp=qb; qp; i++, qp=qp->next);
   return(i);
}
/* procedure 7 */
ATL_tvstrq_t *ATL_ReverseStrQ(ATL_tvstrq_t *oq)
/*
 * RETURNS: Reversed order queue of strings (original is reordered)
 */
{
   ATL_tvstrq_t *nq=NULL, *qp;
   while (oq)
   {
      qp = oq;
      oq = oq->next;
      qp->next = nq;
      nq = qp;
   }
   return(nq);
}

/* procedure 8 */
ATL_tvec_t *ATL_GetTvec(char *name, int N, int nrep, char pre)
{
   ATL_tvec_t *tp;
   int i;

   tp = malloc(sizeof(ATL_tvec_t));
   i = strlen(name);

   tp->name = malloc(sizeof(char)*(strlen(name)+1));
   strcpy(tp->name, name);
   tp->N = N;
   tp->nrep = nrep;
   tp->next = NULL;
   tp->pre = pre;
   i = N;
   if (pre == 's')
      i *= sizeof(char*);
   else if (pre == 'd')
      i *= sizeof(double);
   else
      i *= (pre == 'i') ? sizeof(int) : sizeof(char);
   tp->vp = malloc(N*sizeof(double));
   assert(tp->vp);
   return(tp);
}

/* procedure 9 */
ATL_tvec_t *ATL_KillThisTvec(ATL_tvec_t *tp)
{
   ATL_tvec_t *retp=NULL;

   if (tp)
   {
      if (tp->pre == 's')
      {
         char **sp = tp->vp;
         int i;

         for (i=0; i < tp->N; i++)
            free(sp[i]);
         free(sp);
      }
      else
         free(tp->vp);
      free(tp->name);
      retp = tp->next;
      free(tp);
   }
   return(retp);
}

/* procedure 10 */
void ATL_KillAllTvecs(ATL_tvec_t *tq)
{
   while (tq)
      tq = ATL_KillThisTvec(tq);
}

/* procedure 11 */
void ATL_ReadDoubleTvec(FILE *fpin, int N, double *dp)
{
   int i;

   for (i=0; i < N; i++)
      assert(fscanf(fpin, "%lf\n", dp+i) == 1);
}

/* procedure 12 */
void ATL_ReadIntTvec(FILE *fpin, int N, int *ip)
{
   int i;

   for (i=0; i < N; i++)
      assert(fscanf(fpin, "%d\n", ip+i) == 1);
}

/* procedure 13 */
void ATL_ReadCharTvec(FILE *fpin, int N, char *cp)
{
   int i;

   for (i=0; i < N; i++)
      assert(fscanf(fpin, "%c\n", cp+i) == 1);
}

/* procedure 14 */
void ATL_ReadStringTvec(FILE *fpin, int N, char **sa)
{
   int i;

   for (i=0; i < N; i++)
   {
      char *sp;
      int n, j;
      char ln[512];

      assert(fgets(ln, 512, fpin));
      n = strlen(ln) + 1;
      sa[i] = sp = malloc(sizeof(char)*n);
      assert(sp);
      for (j=0; j < n; j++)
         sp[j] = ln[j];
   }
}

/* procedure 15 */
ATL_tvec_t *ATL_ReadTvec(FILE *fpin)
{
   int N, nrep;
   char nm[64], pre;
   ATL_tvec_t *tp;

   assert(fscanf(fpin, "%s\n", nm) == 1);
   assert(fscanf(fpin, "%d %d %c", &N, &nrep, &pre) == 3);
   tp = ATL_GetTvec(nm, N, nrep, pre);
   if (pre == 'd')
     ATL_ReadDoubleTvec(fpin, N, tp->vp);
   else if (pre == 'i')
     ATL_ReadIntTvec(fpin, N, tp->vp);
   else if (pre == 'c')
     ATL_ReadCharTvec(fpin, N, tp->vp);
   else /* if (pre == 's') */
     ATL_ReadStringTvec(fpin, N, tp->vp);
   return(tp);
}

/* procedure 16 */
ATL_tvec_t *ATL_ReadTvecFile(FILE *fpin, char **cmnt, int *nvec)
/*
 * Reads an entire timing vector file.
 * RETURNS: linked list of timing vectors
 */
{
   int i, n;
   char ln[512];
   ATL_tvec_t *tb, *tp;

   if (!fgets(ln, 512, fpin))
      return(NULL);
   n = strlen(ln);
   while (n > 0 && isspace(ln[n-1]))
      ln[--n] = '\0';;
   *cmnt = malloc(sizeof(char)*n);
   assert(*cmnt);
   strcpy(*cmnt, ln+1);
   assert(fscanf(fpin, " %d\n", nvec) == 1);

   n = *nvec;
   if (n < 1)
      return(NULL);

   tb = tp = ATL_ReadTvec(fpin);
   for (i=1; i < n; i++)
   {
      tp->next = ATL_ReadTvec(fpin);
      tp = tp->next;
   }
   return(tb);
}

/* procedure 17 */
void ATL_WriteTvec(FILE *fpout, ATL_tvec_t *tp)
/*
 * Write a single timing vector to the stream fpout
 */
{
   fprintf(fpout, "%s\n", tp->name);
   fprintf(fpout, "%d %d %c\n", tp->N, tp->nrep, tp->pre);
   if (tp->pre == 'd')
   {
      double *dp = tp->vp;
      const int n = tp->N;
      int i;
      for (i=0; i < n; i++)
         fprintf(fpout, "%le\n", dp[i]);
   }
   else if (tp->pre == 'i')
   {
      int *ip = tp->vp;
      const int n = tp->N;
      int i;
      for (i=0; i < n; i++)
         fprintf(fpout, "%d\n", ip[i]);
   }
   else if (tp->pre == 'c')
   {
      char *cp = tp->vp;
      const int n = tp->N;
      int i;
      for (i=0; i < n; i++)
         fprintf(fpout, "%c\n", cp[i]);
   }
   else /* if (tp->pre == 's') */
   {
      char **sp = tp->vp;
      const int n = tp->N;
      int i;
      for (i=0; i < n; i++)
         fprintf(fpout, "%s\n", sp[i]);
   }
}

/* procedure 18 */
void ATL_WriteTvecs(FILE *fpout, ATL_tvec_t *tp)
/*
 * Writes out a queue of output vectors to a stream that has already had
 * the preample (name, nvec, nrep) written to it
 */
{
   while (tp)
   {
      ATL_WriteTvec(fpout, tp);
      tp = tp->next;
   }
}

/* procedure 19 */
void ATL_WriteTvecFile(FILE *fpout, char *cmnt, int nvec, ATL_tvec_t *tp)
/*
 * Writes the entire output file given a queue of timing vectors
 */
{
   int i;
   fprintf(fpout, "#%s\n", cmnt);
   fprintf(fpout, "%d\n", nvec);
   ATL_WriteTvecs(fpout, tp);
}

/* procedure 20 */
ATL_tvec_t *ATL_FindTvecByName(ATL_tvec_t *tb, char *name)
{
   ATL_tvec_t *tp;
   for (tp=tb; tp && strcmp(tp->name, name); tp = tp->next);
   return(tp);
}

/* procedure 21 */
void ATL_FillCombCHARVecUsingInts
(
   ATL_tvec_t *np,    /* combined vector */
   ATL_tvec_t *ip1,   /* 1st index array (we sort cp1 & cp2 on these ivecs) */
   ATL_tvec_t *ip2,   /* 2nd index array (we sort cp1 & cp2 on these ivecs) */
   ATL_tvec_t *cp1,   /* 1st array to be combined */
   ATL_tvec_t *cp2    /* 2st array to be combined */
)
{
   char *dn = np->vp;
   const char *d1 = cp1->vp, *d2 = cp2->vp;
   const int *s1 = ip1->vp, *s2 = ip2->vp;
   const int n = np->N, n1 = cp1->N, n2 = cp2->N;
   int ic, i1, i2;

   for (ic=i1=i2=0; ic < n; ic++)
   {
      if (i1 < n1)
      {
         if (i2 < n2)  /* both are available for comparison */
         {
            if (s1[i1] <= s2[i2])
               dn[ic] = d1[i1++];
            else
               dn[ic] = d2[i2++];
         }
         else
         {
            assert(i1 < n1);
            dn[ic] = d1[i1++];
         }
      }
      else
      {
         assert(i2 < n2);
         dn[ic] = d2[i2++];
      }
   }
}
/* procedure 22 */
void ATL_FillCombINTVecUsingInts
(
   ATL_tvec_t *np,    /* combined vector */
   ATL_tvec_t *ip1,   /* 1st index array (we sort cp1 & cp2 on these ivecs) */
   ATL_tvec_t *ip2,   /* 2nd index array (we sort cp1 & cp2 on these ivecs) */
   ATL_tvec_t *cp1,   /* 1st array to be combined */
   ATL_tvec_t *cp2    /* 2st array to be combined */
)
{
   int *dn = np->vp;
   const int *d1 = cp1->vp, *d2 = cp2->vp;
   const int *s1 = ip1->vp, *s2 = ip2->vp;
   const int n = np->N, n1 = cp1->N, n2 = cp2->N;
   int ic, i1, i2;

   for (ic=i1=i2=0; ic < n; ic++)
   {
      if (i1 < n1)
      {
         if (i2 < n2)  /* both are available for comparison */
         {
            if (s1[i1] <= s2[i2])
               dn[ic] = d1[i1++];
            else
               dn[ic] = d2[i2++];
         }
         else
         {
            assert(i1 < n1);
            dn[ic] = d1[i1++];
         }
      }
      else
      {
         assert(i2 < n2);
         dn[ic] = d2[i2++];
      }
   }
}
/* procedure 23 */
void ATL_FillCombDOUBLEVecUsingInts
(
   ATL_tvec_t *np,    /* combined vector */
   ATL_tvec_t *ip1,   /* 1st index array (we sort cp1 & cp2 on these ivecs) */
   ATL_tvec_t *ip2,   /* 2nd index array (we sort cp1 & cp2 on these ivecs) */
   ATL_tvec_t *cp1,   /* 1st array to be combined */
   ATL_tvec_t *cp2    /* 2st array to be combined */
)
{
   double *dn = np->vp;
   const double *d1 = cp1->vp, *d2 = cp2->vp;
   const int *s1 = ip1->vp, *s2 = ip2->vp;
   const int n = np->N, n1 = cp1->N, n2 = cp2->N;
   int ic, i1, i2;

   for (ic=i1=i2=0; ic < n; ic++)
   {
      if (i1 < n1)
      {
         if (i2 < n2)  /* both are available for comparison */
         {
            if (s1[i1] <= s2[i2])
               dn[ic] = d1[i1++];
            else
               dn[ic] = d2[i2++];
         }
         else
         {
            assert(i1 < n1);
            dn[ic] = d1[i1++];
         }
      }
      else
      {
         assert(i2 < n2);
         dn[ic] = d2[i2++];
      }
   }
}

/* procedure 24 */
ATL_tvec_t *ATL_CombineTheseVecsUsingInts
(
   ATL_tvec_t *sp1,      /* 1st vector's index array to sort on */
   ATL_tvec_t *sp2,      /* 2nd vector's index array to sort on */
   ATL_tvec_t *cp1,      /* vector to be combined */
   ATL_tvec_t *cp2       /* vector to be combined */
)
{
   ATL_tvec_t *np;
   const pre = cp1->pre;

   assert(sp1->N == cp1->N && sp2->N == cp2->N);
   assert(sp1->pre == sp2->pre && sp1->pre == 'i');
   assert(pre == cp2->pre);
   assert(sp1->nrep == sp2->nrep && sp1->nrep == cp1->nrep &&
          sp1->nrep == cp2->nrep);
   np = ATL_GetTvec(cp1->name, cp1->N+cp2->N, cp1->nrep, pre);
   if (pre == 'i')
      ATL_FillCombINTVecUsingInts(np, sp1, sp2, cp1, cp2);
   else if (pre == 'd')
      ATL_FillCombDOUBLEVecUsingInts(np, sp1, sp2, cp1, cp2);
   else if (pre == 'c')
      ATL_FillCombCHARVecUsingInts(np, sp1, sp2, cp1, cp2);
   return(np);
}

/* procedure 25 */
void ATL_SuffixTvecNames
(
   ATL_tvec_t *tb,   /* list whose names should be suffixed */
   char *suff
)
{
   ATL_tvec_t *tp;
   int isu;

   isu = strlen(suff) + 1;
   for (tp=tb; tp; tp = tp->next)
   {
      int i, inm, n;
      char *sp;

      inm = strlen(tp->name);
      sp = malloc((inm+isu)*sizeof(char));
      for (i=0; i < inm; i++)
         sp[i] = tp->name[i];
      for (n=inm+isu; i < n; i++)
         sp[i] = suff[i-inm];
      free(tp->name);
      tp->name = sp;
   }
}

/* procedure 26 */
ATL_tvec_t *ATL_DupTvec(ATL_tvec_t *t0)
{
   ATL_tvec_t *tp;
   const int N = t0->N;
   int i;

   tp = ATL_GetTvec(t0->name, N, t0->nrep, t0->pre);
   tp->next = NULL;
   if (t0->pre == 'd')
   {
      double *d0 = t0->vp, *d=tp->vp;
      for (i=0; i < N; i++)
         d[i] = d0[i];
   }
   else if (t0->pre == 'i')
   {
      int *d0 = t0->vp, *d=tp->vp;
      for (i=0; i < N; i++)
         d[i] = d0[i];
   }
   else /* pre == 'c' */
   {
      char *d0 = t0->vp, *d=tp->vp;
      assert(t0->pre == 'c');
      for (i=0; i < N; i++)
         d[i] = d0[i];
   }
   return(tp);
}

/* procedure 27 */
void ATL_AdhereTvec
/*
 * Concatonates sp's Tvec entries to the end of dp's.
 */
(
   ATL_tvec_t *dp,      /* destination Tvec */
   ATL_tvec_t *sp       /* source Tvec */
)
{
   ATL_tvec_t *tp;

   if (!sp)
      return;
   assert(dp);
   assert(dp->pre == sp->pre);
   assert(dp->nrep == sp->nrep);
   if (dp->pre == 's')
   {
      const int Nd = dp->N, Ns = sp->N, Nt = Nd+Ns;
      char **sd, **ss;
      int i;
      sd = malloc(Nt * sizeof(char *));
      for (ss=dp->vp, i=0; i < Nd; i++)
         sd[i] = ss[i];
      free(dp->vp);
      dp->vp = sd;
      dp->N = Nt;
      for (ss=sp->vp, i=0; i < Ns; i++)
         sd[i+Nd] = DupString(ss[i]);
   }
   else
   {
      const int Nd = dp->N, Ns = sp->N, Nt = Nd+Ns;
      int sz;
      char *cp;
      if (dp->pre == 'd')
         sz = sizeof(double);
      else
         sz = (dp->pre == 'i') ? sizeof(int) : sizeof(char);
      cp = malloc(Nt*sz);
      assert(cp);
      memcpy(cp, dp->vp, Nd*sz);
      free(dp->vp);
      dp->vp = cp;
      dp->N = Nt;
      memcpy(cp+Nd*sz, sp->vp, Ns*sz);
   }
}

/* procedure 28 */
ATL_tvec_t *ATL_ReverseTvecList  /* RETURNS: reversed list base */
(
   ATL_tvec_t *tb0              /* base of list to reverse */
)
{
   ATL_tvec_t *tb=NULL;
   while (tb0)
   {
      ATL_tvec_t *tp;
      tp = tb0;
      tb0 = tb0->next;
      tp->next = tb;
      tb = tp;
   }
   return(tb);
}

/* procedure 29 */
ATL_tvec_t *ATL_DupNamedVecsFromList  /* returns list of duped vecs */
(
   int N,               /* # of vectors to dup from list to */
   char **names,        /* names of vectors to duplicate */
   ATL_tvec_t *tb0,     /* original list unchanged */
   int order            /* 0: unordered, else: keep in same order as names */
)
{
   ATL_tvec_t *tb=NULL, *tp;
   int i;

   if (!tb0 || N < 1)
      return(NULL);

   for (i=0; i < N; i++)
   {
      tp = ATL_FindTvecByName(tb0, names[i]);
      assert(tp);
      tp = ATL_DupTvec(tp);
      tp->next = tb;
      tb = tp;
   }
   if (order)
      tb = ATL_ReverseTvecList(tb);
   return(tb);
}

/* procedure 30 */
ATL_tvec_t *ATL_RemoveTvecFromList   /* RETURNS: possibly changed base */
(
   ATL_tvec_t *bp,
   ATL_tvec_t *rp
)
{
   ATL_tvec_t *tp = bp, *prev;
   if (!bp)
      return(NULL);
   if (bp == rp)
      return(bp->next);

   prev=tp;
   for (tp=tp->next; tp; tp = tp->next)
   {
      if (tp == rp)
      {
         prev->next = tp->next;
         break;
      }
      prev = tp;
   }
   return(bp);
}

/* procedure 31 */
ATL_tvec_t *ATL_PullNamedVecsFromListWithDups  /* returns list of pulled vecs */
(
   int N,               /* # of vectors to remove from list to */
   char **names,        /* names of vectors to grab */
   ATL_tvec_t **orig    /* original list has names removed */
)
/*
 * Note that the names will be returned in the provided order, with all
 * repeats of a given name contiguous in the list (so in name order, then
 * file order only within names).
 */
{
   ATL_tvec_t *ob=(*orig), *nb=NULL, *tp;
   int i;
   if (!ob || N < 1)
      return(NULL);
   for (i=0; i < N; i++)
   {
/*
 *    Remove all mentions of selected name, add to new list in reverse order
 */
      while ( (tp = ATL_FindTvecByName(ob, names[i])) )
      {
         ob = ATL_RemoveTvecFromList(ob, tp);
         tp->next = nb;
         nb = tp;
      }
   }
   *orig = ob;
   return(ATL_ReverseTvecList(nb));
}

/* procedure 32 */
ATL_tvec_t *ATL_PullNamedVecsFromList  /* returns list of pulled vecs */
(
   int N,               /* # of vectors to remove from list to */
   char **names,        /* names of vectors to grab */
   ATL_tvec_t **orig    /* original list has names removed */
)
/*
 * Note that the names will be returned in the provided order, and that
 * we assume a name only appears once in orig.
 */
{
   ATL_tvec_t *prev, *old=(*orig), *po, *pn, *nb=NULL;
   int i;
   if (!old || N < 1)
      return(NULL);
   for (i=0; i < N; i++)
   {
      prev = NULL;
      po = old;
      while (po)
      {
         if (!strcmp(names[i], po->name))
         {
            if (po == old)
               old = old->next;
            if (nb)
            {
               pn->next = po;
               pn = po;
            }
            else
               pn = nb = po;
            if (prev)
               prev->next = po->next;
            pn->next = NULL;
            break;
         }
         prev = po;
         po = po->next;
      }
   }
   *orig = old;
   return(nb);
}

/* procedure 33 */
int ATL_CountTvecsInList
(
   ATL_tvec_t *tb      /* list to count */
)
{
   int i;
   for (i=0; tb; i++, tb = tb->next);
   return(i);
}

/* procedure 34 */
ATL_tvec_t *ATL_FindLastTvecInList  /* RETURNS: last Tvec in list */
(
   ATL_tvec_t *tb      /* list to look through */
)
{
   if (tb)
   {
      while (tb->next)
         tb = tb->next;
   }
   return(tb);
}

/* procedure 35 */
ATL_tvec_t *ATL_AlphabetizeVecList  /* returns alphabatized list */
(
   int N,              /* # of vectors in list */
   ATL_tvec_t *tb      /* list to alphabetize */
)
/*
 * Alphabatizes N-len tb, and returns new ordered list (old is destroyed).
 */
{
   char **names;
   ATL_tvec_t *tp;
   int i, j;

   names = malloc(sizeof(char*)*N);
   assert(names);
   for (i=0; tp && i < N; i++, tp = tp->next)
      names[i] = tp->name;
   assert(i == N);
/*
 * Sort names using selection sort
 */
   for (j=0; j < N-1; j++)
   {
      for (i=j+1; i < N; i++)
      {
         if (strcmp(names[i], names[j]) < 0)
         {
            char *sp = names[j];
            names[j] = names[i];
            names[i] = sp;
         }
      }
   }
/*
 * Use sorted names to make new alphabetical list, free names and return list
 */
   tb = ATL_PullNamedVecsFromList(N, names, &tb);
   free(names);
   return(tb);
}

/* procedure 36 */
void ATL_CopyStridedVec(char pre, int n, int inc, void *vin, void *vout)
{
   if (pre == 'd')
   {
      int i, j;
      double *in = vin, *out = vout;
      for (j=i=0; i < n; i++, j += inc)
         out[i] = in[j];
   }
   else if (pre == 'i')
   {
      int i, j;
      int *in = vin, *out = vout;
      for (j=i=0; i < n; i++, j += inc)
         out[i] = in[j];
   }
   else if (pre == 'c')
   {
      int i, j;
      char *in = vin, *out = vout;
      for (j=i=0; i < n; i++, j += inc)
         out[i] = in[j];
   }
   else if (pre == 's')
   {
      int i, j;
      char **in = vin, **out = vout;
      for (j=i=0; i < n; i++, j += inc)
      {
         int n;
         n = strlen(in[i]) + 1;
         out[i] = malloc(n*sizeof(char));
         strcpy(out[i], in[j]);
      }
   }
}

/* procedure 37 */
void ATL_PrintTvecElt
(
   FILE *fpout,         /* stream to print to */
   ATL_tvec_t *tp,      /* vector to print from */
   int idx              /* index in vector to print */
)
{
   if (tp->pre == 'd')
   {
      double *p = tp->vp;
      fprintf(fpout, "%e", p[idx]);
   }
   else if (tp->pre == 'i')
   {
      int *p = tp->vp;
      fprintf(fpout, "%12d", p[idx]);
   }
   else if (tp->pre == 'c')
   {
      char *p = tp->vp;
      fprintf(fpout, "%c", p[idx]);
   }
   else if (tp->pre == 's')
   {
      char **p = tp->vp;
      fprintf(fpout, "%s", p[idx]);
   }
}

/* procedure 38 */
void ATL_PrintTvecsInRow
(
   FILE *fpout,         /* stream to print to */
   ATL_tvec_t *tb,      /* list of vectors to print (rowwise) */
   char *sep,           /* string to print between elements */
   int head             /* should we print headers? */
)
{
   ATL_tvec_t *tp;
   const int N = tb->N;
   int i;

   if (!tb)
      return;
   for (tp=tb->next; tp; tp = tp->next)
      assert(tp->N >= N);

   if (head)
      for (tp=tb; tp; tp = tp->next)      /* loop over columns  */
         fprintf(fpout, "%12.12s%s", tp->name, sep);
   fprintf(fpout, "%s\n", sep);
   for (i=0; i < N; i++)                        /* loop over rows of vectors */
   {
      for (tp=tb; tp->next; tp = tp->next)      /* loop over columns  */
      {
         ATL_PrintTvecElt(fpout, tp, i);
         fprintf(fpout, "%s", sep);
      }
      ATL_PrintTvecElt(fpout, tp, i);
      fprintf(fpout, "%s\n", sep);
   }
}

/* procedure 39 */
ATL_tvec_t *ATL_GetRep1Tvec
(
   ATL_tvec_t *tin,    /* vector to split */
   int istart          /* which repetition to grab */
)
{
   char *name;
   ATL_tvec_t *t1=NULL;
   int n;
   const int Nr = tin->N / tin->nrep;
   const char pre = tin->pre;

   assert(tin->nrep < 10000000);
   n = strlen(tin->name) + 9;
   name = malloc(n*sizeof(char));
   assert(name);
   sprintf(name, "%s_%d", tin->name, istart);
   t1 = ATL_GetTvec(name, Nr, 1, pre);
   free(name);
   if (pre == 'd')
      ATL_CopyStridedVec(pre, Nr, tin->nrep,((double*)(tin->vp))+istart,t1->vp);
   else if (pre == 'i')
      ATL_CopyStridedVec(pre, Nr, tin->nrep, ((int*)(tin->vp))+istart, t1->vp);
   else if (pre == 'c')
      ATL_CopyStridedVec(pre, Nr, tin->nrep, ((char*)(tin->vp))+istart, t1->vp);
   else if (pre == 's')
      ATL_CopyStridedVec(pre, Nr, tin->nrep, ((char**)(tin->vp))+istart,t1->vp);

   return(t1);
}

/* procedure 40 */
ATL_tvec_t *ATL_SplitRepsTvec  /* returns Q of sep vecs for each rep */
(
   ATL_tvec_t *tin     /* vector to split */
)
{
   int i;
   const int nrep = tin->nrep;
   ATL_tvec_t *tb, *tp;

   tp = tb = ATL_GetRep1Tvec(tin, 0);
   for (i=1; i < nrep; i++)
   {
      tp->next = ATL_GetRep1Tvec(tin, i);
      tp = tp->next;
   }
   return(tb);
}

/* procedure 41 */
ATL_tvec_t *ATL_GetStatTvecsDOUBLE  /* returns Q of stat vectors: <,+,> */
(
   ATL_tvec_t *tin     /* vector to get stats for */
)
{
   char *name;
   int i, j, n, N;
   ATL_tvec_t *tavg, *tmin, *tmax;
   double *dmin, *dmax, *davg;

   n = strlen(tin->name) + 5;
   name = malloc(sizeof(char)*n);
   assert(name);
   N = tin->N / tin->nrep;

   sprintf(name, "%s_min", tin->name);
   tmin = ATL_GetTvec(name, N, 1, tin->pre);
   sprintf(name, "%s_avg", tin->name);
   tavg = ATL_GetTvec(name, N, 1, tin->pre);
   sprintf(name, "%s_max", tin->name);
   tmax = ATL_GetTvec(name, N, 1, tin->pre);
   free(name);
   dmin = tmin->vp;
   dmax = tmax->vp;
   davg = tavg->vp;

   n = tin->nrep;
   for (j=0; j < N; j++)
   {
      double min, max, sum, *din;

      din = ((double*)(tin->vp)) + j*n;
      min = max = sum = *din;
      for (i=1; i < n; i++)
      {
         if (din[i] < min)
            min = din[i];
         if (din[i] > max)
            max = din[i];
         sum += din[i];
      }
      dmin[j] = min;
      dmax[j] = max;
      davg[j] = sum / n;
   }
   tmin->next = tavg;
   tavg->next = tmax;
   return(tmin);
}
/*
*******************************************************************************
*                   OVERVIEW OF FUNCTIONS IN THIS FILE                        *
*******************************************************************************
* ATL_tvec_t *ATL_GetStatTvecsDOUBLE(ATL_tvec_t *tin)
*    returns list of stat vectors: <,+,>
* ATL_tvec_t *ATL_SplitRepsTvec(ATL_tvec_t *tin)
*    Returns list of seperate vecs for each rep (tin unchanged)
* ATL_tvec_t *ATL_GetRep1Tvec(ATL_tvec_t *tin, int istart)
*    Return new nrep=1 tvec wt only repitition istart represented
* void ATL_PrintTvecsInRow(FILE *fpout, ATL_tvec_t *tb, char *sep, int head)
*    Prints tb, 1 elt per line tvec cols seperated by sep
* void ATL_PrintTvecElt(FILE *fpout, ATL_tvec_t *tp, int idx)
*    Print element idx of tp->vp
* void ATL_CopyStridedVec(char pre, int n, int inc, void *vin, void *vout)
*    copies every inc elts of vin to vout
* ATL_tvec_t *ATL_AlphabetizeVecList(int N, ATL_tvec_t *tb)
*    Alphabatizes N-len tb, and returns new ordered list (old destroyed)
* ATL_tvec_t *ATL_FindLastTvecInList(ATL_tvec_t *tb)
*    Returns last tvec in list
* int ATL_CountTvecsInList(ATL_tvec_t *tb)
*    Returns number of tvecs in tb
* ATL_tvec_t *ATL_PullNamedVecsFromList(int N, char **names, ATL_tvec_t bp)
*    Returns list of first instance of named tvecs, which are removed from bp
* ATL_tvec_t *ATL_PullNamedVecsFromListWithDups(int N, char **names, tvec_t bp)
*    Returns new list of all named tvecs, which are removed from bp
* ATL_tvec_t *ATL_RemoveTvecFromList(ATL_tvec_t *bp, ATL_tvec_t *rp)
*    Returns new bp after removing rp from bp
* ATL_tvec_t *ATL_DupNamedVecsFromList(int N, char** nm, tvec_t *tb, int order)
*    Returns dups of named vectors from tb; repeat names not duped!
* ATL_tvec_t *ATL_ReverseTvecList(ATL_tvec_t *tb)
*    Reverse order in list tb
* void ATL_AdhereTvec(ATL_tvec_t *dp, ATL_tvec_t *sp)
*    Stick data in sp->vp to end of dp->vp, sp not changed
* ATL_tvec_t *ATL_DupTvec(ATL_tvec_t *t0)
*    Return newly allocated tvec identical to t0, except ->next is NULL
* void ATL_SuffixTvecNames(ATL_tvec_t *tb, char *suff)
*    suffix all names in tb with suff
* ATL_tvec_t *ATL_CombineTheseVecsUsingInts
*    Merges two tvecs, go read the function
* void ATL_FillCombDOUBLEVecUsingInts
*    Merges two tvecs, go read the function
* void ATL_FillCombINTVecUsingInts
*    Merges two tvecs, go read the function
* void ATL_FillCombCHARVecUsingInts
*    Merges two tvecs, go read the function
* ATL_tvec_t *ATL_FindTvecByName(ATL_tvec_t *tb, char *name)
*    Returns first tvec with name name (NULL if not found)
* void ATL_WriteTvecFile(FILE *fpout, char *cmnt, int nvec, ATL_tvec_t *tp)
*    Writes the entire output file given a queue of timing vectors
* void ATL_WriteTvecs(FILE *fpout, ATL_tvec_t *tp)
*    Writes list of tvecs assuming stream already has preamble in it
* void ATL_WriteTvec(FILE *fpout, ATL_tvec_t *tp)
*    Write a single timing vector to the stream fpout
* ATL_tvec_t *ATL_ReadTvecFile(FILE *fpin, char **cmnt, int *nvec)
*    Return all tvecs from file fpin (fpin at beginning of file!)
* ATL_tvec_t *ATL_ReadTvec(FILE *fpin)
*    Returns one tvec read in from fpin
* ATL_ReadStringTvec(FILE *fpin, int N, char **sa)
*    Reads in N strings from fpin into sa
* ATL_ReadCharTvec(FILE *fpin, int N, char *cp)
*    Reads in N chars from fpin into cp
* void ATL_ReadIntTvec(FILE *fpin, int N, int *ip)
*    Reads in N int from fpin into ip
* void ATL_ReadDoubleTvec(FILE *fpin, int N, double *dp)
*    Reads in N doubles from fpin into dp
* void ATL_KillAllTvecs(ATL_tvec_t *tq)
*    Frees all Tvecs (and internals) in list
* ATL_tvec_t *ATL_KillThisTvec(ATL_tvec_t *tp)
*    Frees all internal ptrs, frees arg tvec (tp), returns tp->next
* ATL_tvec_t *ATL_GetTvec(char *name, int N, int nrep, char pre)
*    Returns allocated Tvec struct, next is NULL, vp alloc, but not init
* ATL_tvstrq_t *ATL_ReverseStrQ(ATL_tvstrq_t *oq)
*    Returns oq in reverse order
* int ATL_CountStrNodes(ATL_tvstrq_t *qb)
*    Returns number of nodes in qb
* char **ATL_StrQ2Arr(ATL_tvstrq_t *qb, int *N)
*    Kills qb to make N+1-length NULL-term array of strings instead
* ATL_tvstrq_t *RemoveStrNodeFromQ(ATL_tvstrq_t *b, ATL_tvstrq_t *p)
*    removes p from Q b, sets p->next=NULL, returns new b
* void ATL_KillStrQ(ATL_tvstrq_t *kb)
*    Deallocates entire string queue
* ATL_tvstrq_t *ATL_KillStrNode(ATL_tvstrq_t *kp)
*    Deallocates node & string, returns next
* ATL_tvstrq_t *ATL_GetStrNode(char *name)
*    Allocates node & new string dup, returns node ptr
*******************************************************************************
*******************************************************************************
*/

#endif  /* end ifdef multi-inclusion guard */
