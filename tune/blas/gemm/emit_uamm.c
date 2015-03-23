/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 2013 R. Clint Whaley
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
#define ATL_JKMDEF 1
#include "atlas_mmparse.h"
#include "atlas_sys.h"
   static int UID=0, UIL=1;
void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set type/precision prefix (d) \n");
   fprintf(stderr, "   -d <outdir>: directory to dump files to\n");
   fprintf(stderr, "   -i <infile> : can be repeated for multiple files\n");
   fprintf(stderr, "   -k <unique K cleanup index file> : \n");
   fprintf(stderr, "   -K <K cleanup by NB file> \n");
   fprintf(stderr, "   -r <rank-K kernel file> \n");
   fprintf(stderr, "   -s <square-case kernel file>\n");

   fprintf(stderr,
      "   -I <ID> : unique non-negative ID for header/kern files\n");
   exit(ierr ? ierr : -1);
}

ATL_mmnode_t *GetFlags(int nargs, char **args, char *PRE, char **DOUT,
                       char **UKIN, char **KCIN, char **RKIN, char **SQIN)
{
   int i, j=0, n, k;
   char pre='d';
   *SQIN = *RKIN = *UKIN = *KCIN = *DOUT = NULL;
   ATL_mmnode_t *mmb=NULL, *mmp, *mp;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'p':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        pre = tolower(args[i][0]);
        assert(pre == 's' || pre == 'd' || pre == 'z' || pre == 'c');
        break;
      case 's':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *SQIN = DupString(args[i]);
        break;
      case 'k':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *UKIN = DupString(args[i]);
        break;
      case 'K':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *KCIN = DupString(args[i]);
        break;
      case 'I':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        UID = atol(args[i]);
        for (k=10; k <= UID; k *= 10)
           UIL++;
        break;
      case 'd':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *DOUT = DupString(args[i]);
        break;
      case 'i':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        mmp = ReadMMFile(args[i]);
        if (mmb)
        {
           ATL_mmnode_t *mp;
           for (mp=mmb; mp->next; mp = mp->next);
           mp->next = mmp;
        }
        else
           mmb = mmp;
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   *PRE = pre;
   if (!(*DOUT))
   {
      *DOUT = DupString("dMake_amm");
      (*DOUT)[0] = pre;
   }
   return(mmb);
}

char *GetVecStr(char pre, int vlen)
{
   if (vlen == 1)
      return("scalar");
   #ifdef ATL_AVX
      if (pre == 'd' || pre == 'z')
      {
         if (vlen == 4)
            return("avx");
         else if (vlen == 2)
            return("sse");
      }
      else if (pre == 's' || pre == 'c')
      {
         if (vlen == 8)
            return("avx");
         else if (vlen == 4)
            return("sse");
      }
   #elif defined(ATL_SSE1)
      #ifdef ATL_SSE2
         if ((pre == 'd' || pre == 'z') && vlen == 2)
               return("sse");
      #endif
      if ((pre == 's' || pre == 'c') && vlen == 4)
         return("sse");
   #endif
/*
 * Any vector length > 1 that isn't one of our known cases uses gnuvec
 */
   return("gvec");
}

void FillInGenStrings
(
   char pre,
   ATL_mmnode_t *mmb,  /* queue to look through */
   char *dir           /* output directory to generate into */
)
/*
 * Creates GenString for any ID=0 in mmb
 */
{
   ATL_mmnode_t *mp;
   char *frm="make gen_amm_%s lat=%d mu=%d nu=%d ku=%d kb=%d vlen=%d rt=%s/%s kmaj=%d";
   int frml;

   frml = strlen(frm);
   for (mp=mmb; mp; mp = mp->next)
   {
      if (mp->ID == 0)  /* is generated file */
      {
         int ll;
         char *sp, *vec;
         int vlen = mp->vlen;

         assert(mp->rout);  /* should have been filled in by search */
         ll = frml + strlen(dir) + strlen(mp->rout) + 8;
         sp = malloc(ll*sizeof(char));
         assert(sp);

         vec = GetVecStr(pre, mp->vlen);

         sprintf(sp, frm, vec, mp->lat, mp->mu, mp->nu,
                 mp->ku, FLAG_IS_SET(mp->flag, MMF_KRUNTIME) ? 0 : mp->kbB,
                 vlen, dir, mp->rout, mp->kmaj);
         mp->genstr = sp;
      }
   }
}

void PrintBegBlock(char pre, ATL_mmnode_t *mmb, char *nam, FILE *fp)
{
   ATL_mmnode_t *mp;
   char PRE = toupper(pre);
   int i;

   if (nam)
   {
      fprintf(fp, "#ifndef ATLAS_%cUAMM_%s_H\n   #define ATLAS_%cUAMM_%s_H\n\n",
              PRE, nam, PRE, nam);
      fprintf(fp, "#include \"atlas_amm.h\"\n");
   }
   else
      fprintf(fp, "#ifndef ATLAS_%cUAMM_H\n   #define ATLAS_%cUAMM_H\n\n",
              PRE, PRE);
/*
 * Count mmb, and print def of NCASES
 */
   if (!nam || strstr(nam, "RANKK") == NULL)
   {
      for (mp=mmb,i=0; mp; i++, mp = mp->next);

      fprintf(fp, "#ifdef ATL_UAMM_NCASES\n");
      fprintf(fp, "   #if ATL_UAMM_NCASES != %d\n", i);
      fprintf(fp, "      #error \"NCASES MISMATCH!\"\n");
      fprintf(fp, "   #endif\n");
      fprintf(fp, "#else\n");
      fprintf(fp, "   #define ATL_UAMM_NCASES %d\n", i);
      fprintf(fp, "#endif\n");
   }
}

char *GetHName(char pre, char *outd, char *bnam)
{
   int i, NOBASE=0;
   char *fnam;
   if (!bnam)
   {
      NOBASE = 1;
      bnam = "";
   }
   i = strlen(outd) + strlen(bnam) + 16+UIL;

   fnam = malloc(i*sizeof(char));
   assert(fnam);
   if (NOBASE)
      sprintf(fnam, "%s/atlas_%cu%damm.h", outd, pre, UID);
   else
      sprintf(fnam, "%s/atlas_%cu%damm_%s.h", outd, pre, UID, bnam);
   return(fnam);
}

FILE *StandHStart(char pre, ATL_mmnode_t *mmb, char *outd, char *bnam)
{
   char *fnam;
   FILE *fp;
   int i;

   assert(outd);
   fnam = GetHName(pre, outd, bnam);
   fp = fopen(fnam, "w");
   assert(fp);
   if (bnam)
   {
      for (i=0; bnam[i]; i++)
         fnam[i] = toupper(bnam[i]);
      fnam[i] = '\0';
      PrintBegBlock(pre, mmb, fnam, fp);
   }
   else
      PrintBegBlock(pre, mmb, NULL, fp);
   free(fnam);
   return(fp);
}

static int Mylcm(const int M, const int N)
/*
 * Returns least common multiple (LCM) of two positive integers M & N by
 * computing greatest common divisor (GCD) and using the property that
 * M*N = GCD*LCM.
 */
{
   register int tmp, max, min, gcd=0;

   if (M != N)
   {
      if (M > N) { max = M; min = N; }
      else { max = N; min = M; }
      if (min > 0)  /* undefined for negative numbers */
      {
         do  /* while (min) */
         {
            if ( !(min & 1) ) /* min is even */
            {
               if ( !(max & 1) ) /* max is also even */
               {
                  do
                  {
                     min >>= 1;
                     max >>= 1;
                     gcd++;
                     if (min & 1) goto MinIsOdd;
                  }
                  while ( !(max & 1) );
               }
               do min >>=1 ; while ( !(min & 1) );
            }
/*
 *          Once min is odd, halve max until it too is odd.  Then, use
 *          property that gcd(max, min) = gcd(max, (max-min)/2)
 *          for odd max & min
 */
MinIsOdd:
            if (min != 1)
            {
               do  /* while (max >= min */
               {
                  max -= (max & 1) ? min : 0;
                  max >>= 1;
               }
               while (max >= min);
            }
            else return( (M*N) / (1<<gcd) );
            tmp = max;
            max = min;
            min = tmp;
         }
         while(tmp);
      }
      return( (M*N) / (max<<gcd) );
   }
   else return(M);
}

void GenAmmSum(char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *rkb, char *outd)
{
   ATL_mmnode_t *mp, *p66;
   char *fnam;
   FILE *fp;
   char *type = "unsigned short";
   int i, n, maxb, maxNB, maxMB, maxKB, maxkmaj;
   char PRE = toupper(pre), bc[3] = {'M', 'N', 'K'};
   double mfB;

   fp = StandHStart(pre, mmb, outd, "sum");
   maxkmaj = maxNB = maxKB = maxMB = 0;
   mfB = 0.0;
   for (n=0,mp=mmb; mp; n++, mp = mp->next)
   {
      maxkmaj = Mmax(maxkmaj, mp->kmaj);
      maxMB = Mmax(maxMB, mp->mbB);
      maxNB = Mmax(maxNB, mp->nbB);
      maxKB = Mmax(maxKB, mp->kbB);
      mfB = Mmax(mfB, mp->mflop[0]);
   }
   if (maxkmaj == 1)
      maxkmaj = 0;
   maxb = Mmax(maxMB, maxNB);
   maxb = Mmax(maxb, maxKB);
   fprintf(fp, "\n#define ATL_UAMM_MAXMB %d\n", maxMB);
   fprintf(fp, "#define ATL_UAMM_MAXNB %d\n", maxNB);
   fprintf(fp, "#define ATL_UAMM_MAXKB %d\n", maxKB);
   fprintf(fp, "#define ATL_UAMM_MAXKMAJ %d\n\n", maxkmaj);

   for (mp=mmb; mp && mp->next; mp = mp->next);
   assert(mp);
   fprintf(fp, "#define ATL_AMM_LMU %d\n", mp->mu);
   fprintf(fp, "#define ATL_AMM_LNU %d\n", mp->nu);
   fprintf(fp, "#define ATL_AMM_LKU %d\n", mp->ku);
   fprintf(fp, "#define ATL_AMM_LLCMMN %d\n\n", Mylcm(mp->mu, mp->nu));
   fprintf(fp, "#define ATL_AMM_LLCMU %d\n\n",
           Mylcm(Mylcm(mp->mu, mp->nu),mp->ku));
/*
 * Find smallest case achieving 2/3 of maximal performance
 */
   for (i=0,mp=mmb; mp && mp->mflop[0]*1.5 < mfB; i++, mp = mp->next);
   assert(mp);
   fprintf(fp, "#define ATL_UAMM_66IDX %d\n", i);
   fprintf(fp, "#define ATL_UAMM_66MB %d\n", mp->mbB);
   fprintf(fp, "#define ATL_UAMM_66NB %d\n", mp->nbB);
   fprintf(fp, "#define ATL_UAMM_66KB %d\n", mp->kbB);
   fprintf(fp, "#define ATL_AMM_66LCMMN %d\n\n", Mylcm(mp->mu, mp->nu));
   fprintf(fp, "#define ATL_AMM_66LCMU %d\n\n",
           Mylcm(Mylcm(mp->mu, mp->nu),mp->ku));
   fprintf(fp, "#define ATL_UAMM_66RATIO %1.4lf\n\n", mp->mflop[0]/mfB);
/*
 * Find smallest case achieving 98% of maximal performance
 */
   for (i=0,mp=mmb; mp && mp->mflop[0] < 0.98*mfB; i++, mp = mp->next);
   assert(mp);
   fprintf(fp, "#define ATL_UAMM_98IDX %d\n", i);
   fprintf(fp, "#define ATL_UAMM_98MB %d\n", mp->mbB);
   fprintf(fp, "#define ATL_UAMM_98NB %d\n", mp->nbB);
   fprintf(fp, "#define ATL_UAMM_98KB %d\n", mp->kbB);
   fprintf(fp, "#define ATL_AMM_98LCMMN %d\n\n", Mylcm(mp->mu, mp->nu));
   fprintf(fp, "#define ATL_AMM_98LCMU %d\n\n",
           Mylcm(Mylcm(mp->mu, mp->nu),mp->ku));
   fprintf(fp, "#define ATL_UAMM_98RATIO %1.4lf\n\n", mp->mflop[0]/mfB);
   assert(rkb == NULL);

   fprintf(fp, "#define ATL_AMMFLG_KRUNTIME(flg_) ((flg_) & 1)\n");
   fprintf(fp, "#define ATL_AMMFLG_KMAJOR(flg_) ((flg_) & 2)\n");

   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenPerfFile(char pre, ATL_mmnode_t *mmb, char *outd, char *nm)
{
   ATL_mmnode_t *mp;
   #define NTHRSH 11
   int THRSH[NTHRSH] = {25, 33, 50, 66, 75, 80, 85, 90, 95, 98, 99};
   int idxT[NTHRSH];
   ATL_mmnode_t *mpT[NTHRSH];
   char *fnam;
   FILE *fp;
   char *type = "float";
   double mfMax=0.0;
   int i, j, n, maxb, maxNB, maxMB, maxKB, maxkmaj, idxMax=0;
   char PRE = toupper(pre), bc[3] = {'M', 'N', 'K'};

   for (i=0; i < NTHRSH; i++)
      mpT[i] = NULL;
   fp = StandHStart(pre, mmb, outd, nm);
   for (n=0,mp=mmb; mp; n++, mp = mp->next)
   {
      if (mp->mflop[0] > mfMax)
      {
         mfMax = mp->mflop[0];
         idxMax = n;
      }
   }
   fprintf(fp, "#define ATL_UAMM_MAXMFLOP %le /* (%.2f)*/ \n",
           mfMax, mfMax);
   fprintf(fp, "#define ATL_UAMM_MAXMFLOPIDX %d\n\n", idxMax);
   for (n=0,mp=mmb; mp; mp = mp->next, n++)
   {
      double mf = mp->mflop[0] / mfMax;
      for (i=0; i < NTHRSH; i++)
      {
         if (!mpT[i] && THRSH[i]*0.01*mfMax < mp->mflop[0])
         {
            mpT[i] = mp;
            idxT[i] = n;
         }
      }
   }
   for (i=0; i < NTHRSH; i++)
   {
      mp = mpT[i];
      fprintf(fp, "#define ATL_UAMM_%dLCMU %d\n", THRSH[i],
              Mylcm(Mylcm(mp->mu,mp->nu),mp->ku));
      fprintf(fp, "#define ATL_UAMM_%dLCMMN %d\n", THRSH[i],
              Mylcm(mp->mu,mp->nu));
      fprintf(fp, "#define ATL_UAMM_%dKB %d\n", THRSH[i],
              mp->kbB);
      fprintf(fp, "#define ATL_UAMM_%dNB %d\n", THRSH[i],
              mp->nbB);
      fprintf(fp, "#define ATL_UAMM_%dMB %d\n", THRSH[i],
              mp->mbB);
      fprintf(fp, "#define ATL_UAMM_%dIDX %d\n", THRSH[i],
              idxT[i]);
   }
   fprintf(fp, "\n");

   fprintf(fp, "static const float ATL_UAMM_PERF[%d] =", n);
   fprintf(fp, "   /* %% of performance of best kernel */\n{\n");
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
      fprintf(fp, "   %f%c  /* IDX=%d, KB=%d */\n", mp->mflop[0]/mfMax,
              (mp->next)?',':' ', j, mp->kbB);
   fprintf(fp, "};\n\n");

   fprintf(fp, "static const float ATL_UAMM_SPDUPNXT[%d] =", n);
   fprintf(fp, "   /* speedup of next higher NB */\n{\n");
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      double mf = (mp->next) ? mp->next->mflop[0] : mp->mflop[0];
      fprintf(fp, "   %f%c  /* IDX=%d, KB=%d vs. %d */\n", mf/mp->mflop[0],
              (mp->next)?',':' ', j, mp->kbB, mp->next?mp->next->kbB:mp->kbB);
   }
   fprintf(fp, "};\n\n");

   fprintf(fp, "#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenBlockingFile(char pre, ATL_mmnode_t *mmb, char *outd, char *nm)
{
   ATL_mmnode_t *mp;
   char *fnam;
   FILE *fp;
   char *type = "unsigned short";
   int i, n, maxb, maxNB, maxMB, maxKB, maxkmaj;
   char PRE = toupper(pre), bc[3] = {'M', 'N', 'K'};

   fp = StandHStart(pre, mmb, outd, nm);
   maxkmaj = maxNB = maxKB = maxMB = 0;
   for (n=0,mp=mmb; mp; n++, mp = mp->next)
   {
      maxkmaj = Mmax(maxkmaj, mp->kmaj);
      maxMB = Mmax(maxMB, mp->mbB);
      maxNB = Mmax(maxNB, mp->nbB);
      maxKB = Mmax(maxKB, mp->kbB);
   }
   if (maxkmaj == 1)
      maxkmaj = 0;
   maxb = Mmax(maxMB, maxNB);
   maxb = Mmax(maxb, maxKB);
   fprintf(fp, "#define ATL_UAMM_MAXMB %d\n", maxMB);
   fprintf(fp, "#define ATL_UAMM_MAXNB %d\n", maxNB);
   fprintf(fp, "#define ATL_UAMM_MAXKB %d\n", maxKB);
   fprintf(fp, "#define ATL_UAMM_MAXKMAJ %d\n", maxkmaj);
   fprintf(fp, "\n");

   if (maxb <= 255)
      type = "unsigned char";
   for (i=0; i < 3; i++)
   {
      int j;
      fprintf(fp, "static const %s ATL_UAMM_%cBs[%d] =\n{\n",
              type, bc[i], n);
      for (j=0,mp=mmb; mp; j++,mp = mp->next)
      {
         int b;
         if (bc[i] == 'M')
            b = mp->mbB;
         else if (bc[i] == 'N')
            b = mp->nbB;
         else if (bc[i] == 'K')
            b = mp->kbB;
         if (mp->next)
            fprintf(fp, "%8d,  /* index %d */\n", b, j);
         else
            fprintf(fp, "%8d   /* index %d */\n", b, j);
      }
      fprintf(fp, "};\n\n");
   }
   for (i=0; i < 3; i++)
   {
      int j;
      fprintf(fp, "static const %s ATL_UAMM_%cUs[%d] =\n{\n",
              type, bc[i], n);
      for (j=0,mp=mmb; mp; j++,mp = mp->next)
      {
         int b;
         if (bc[i] == 'M')
            b = mp->mu;
         else if (bc[i] == 'N')
            b = mp->nu;
         else if (bc[i] == 'K')
         {
            if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
               b = mp->ku;
            else
               b = (mp->kmaj > 1) ? mp->ku : mp->kbB;
         }
         if (mp->next)
            fprintf(fp, "%8d,  /* index %d */\n", b, j);
         else
            fprintf(fp, "%8d   /* index %d */\n", b, j);
      }
      fprintf(fp, "};\n\n");
   }
   fprintf(fp, "static const %s ATL_UAMM_KBMINs[%d] =\n{\n", type, n);
   for (i=0,mp=mmb; mp; i++,mp = mp->next)
   {
      if (mp->next)
         fprintf(fp, "%8d,  /* index %d */\n", mp->kbmin, i);
      else
         fprintf(fp, "%8d   /* index %d */\n", mp->kbmin, i);
   }
   fprintf(fp, "};\n\n");
   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenFlagH(char pre, ATL_mmnode_t *mmb, char *outd, char *nm)
{
   FILE *fp;
   int j, n;
   ATL_mmnode_t *mp;

   fp = StandHStart(pre, mmb, outd, nm);

   for (n=0,mp=mmb; mp; n++,mp = mp->next);

   fprintf(fp, "static const unsigned char ATL_UAMM_KFLAG[%d] =\n{\n", n);
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      unsigned char flag=FLAG_IS_SET(mp->flag, MMF_KRUNTIME) ? 1 : 0;
      if (mp->kmaj)
         flag |= 2;
      if (mp->next)
         fprintf(fp, "%6d,  /* index %d */\n", flag, j);
      else
         fprintf(fp, "%6d   /* index %d */\n", flag, j);
   }
   fprintf(fp, "};\n\n");
   fprintf(fp, "#define ATL_AMM_KRUNTIME(idx_) (ATL_AMM_KFLAG[idx_] & 1)\n");
   fprintf(fp, "#define ATL_AMM_KMAJOR(idx_) (ATL_AMM_KFLAG[idx_] & 2)\n");
   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void SpewForthC2MProto(char pre, FILE *fp0, FILE *fp1, int mu, int nu)
{
   char ac[3] = {'1', 'n', 'X'};
   char bc[4] = {'0', '1', 'n', 'X'};
   int ia, ib;
   for (ia=0; ia < 3; ia ++)
   {
      for (ib=0; ib < 4; ib++)
      {
         fprintf(fp0, "void ATL_%cu%dablk2cmat_%dx%d_a%c_b%c\n",
                 pre, UID, mu, nu, ac[ia], bc[ib]);
         if (pre == 'z' || pre == 'c')
            fprintf(fp0, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const TYPE*,const SCALAR,TYPE *,ATL_CSZT);\n");
         else
            fprintf(fp0, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const SCALAR,TYPE *,ATL_CSZT);\n");
         fprintf(fp1, "void ATL_%cu%dcmat2ablk_%dx%d_a%c_b%c\n",
                 pre, UID, mu, nu, ac[ia], bc[ib]);
         if (pre == 'z' || pre == 'c')
            fprintf(fp1, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,TYPE*);\n");
         else
            fprintf(fp1, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,const SCALAR,TYPE*);\n");
      }
   }
}

void SpewForthC2BDecl(char pre, ATL_mmnode_t *mmb, FILE *fp, char *rt,
                      char alp, char bet)
{
   ATL_mmnode_t *mp;
   int j;

   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      fprintf(fp, "   ATL_%cu%d%s_%dx%d_a%c_b%c",
              pre, UID, rt, mp->mu,mp->nu, alp, bet);
      if (mp->next)
         fprintf(fp, ",  /* index %d */\n", j);
      else
         fprintf(fp, "   /* index %d */\n", j);
      }
      fprintf(fp, "};\n\n");
}

void GenC2BLK(char pre, ATL_mmnode_t *mmb, char *outd, char *suff)
{
   FILE *fp0, *fp1;
   ATL_mmnode_t *mp;
   int ia, ib;
   char ac[3] = {'1', 'n', 'X'};
   char bc[4] = {'0', '1', 'n', 'X'};
   char *fnam;

   if (!suff)
   {
      fp0 = StandHStart(pre, mmb, outd, "ablk2cmat");
      fp1 = StandHStart(pre, mmb, outd, "cmat2ablk");
   }
   else
   {
      fnam = malloc(sizeof(char)*(strlen(suff) + 10));
      assert(fnam);
      strncpy(fnam, "ablk2cmat", 9);
      strcpy(fnam+9, suff);
      fp0 = StandHStart(pre, mmb, outd, fnam);
      strncpy(fnam, "cmat2ablk", 9);
      fp1 = StandHStart(pre, mmb, outd, fnam);
   }
   fprintf(fp0, "\n");
   fprintf(fp1, "\n");
/*
 * Crank out prototypes
 */
   SpewForthC2MProto(pre, fp0, fp1, mmb->mu, mmb->nu);
   for (mp=mmb->next; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      const int mu=mp->mu, nu=mp->nu;
      for (p=mmb; p->mu != mu || p->nu != nu; p = p->next);
      if (p == mp)  /* first occurance of this mu,nu pair */
         SpewForthC2MProto(pre, fp0, fp1, mp->mu, mp->nu);
   }
   fprintf(fp0, "\n");
   fprintf(fp1, "\n");
/*
 * Now, crank out funcptr arrays
 */
   for (ia=0; ia < 3; ia ++)
   {
      for (ib=0; ib < 4; ib++)
      {
         fprintf(fp0,
            "static const ablk2cmat_t ATL_UAMM_BLK2C_a%c_b%c[ATL_UAMM_NCASES] =\n{\n",
                 ac[ia], bc[ib]);
         SpewForthC2BDecl(pre, mmb, fp0, "ablk2cmat", ac[ia], bc[ib]);
         fprintf(fp1,
            "static const cmat2ablk_t ATL_UAMM_C2BLK_a%c_b%c[ATL_UAMM_NCASES] =\n{\n",
                 ac[ia], bc[ib]);
         SpewForthC2BDecl(pre, mmb, fp1, "cmat2ablk", ac[ia], bc[ib]);
      }
   }
   fprintf(fp0, "\n#endif  /* end include file guard */\n");
   fclose(fp0);
   fprintf(fp1, "\n#endif  /* end include file guard */\n");
   fclose(fp1);
}

void SpewForthRevCpProto(char pre, FILE *fp, char alp, int u, int kmaj)
{
   const int G = (pre == 'c' || pre == 'z') ? 2 : 1;
   const char *cst[2] = {"", "C"};
   int g;

   for (g=0; g < G; g++)
   {
      if (kmaj > 1)
         fprintf(fp, "void ATL_%cu%dam2cm_a%c_%dx%d%s\n",
                 pre, UID, alp, kmaj, u, cst[g]);
      else
         fprintf(fp, "void ATL_%cu%dam2cm_a%c_%d%s\n",pre, UID, alp, u, cst[g]);
      if (pre == 'z' || pre == 'c')
         fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*,const TYPE*);\n");
      else
         fprintf(fp,
         "   (ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*);\n");
      if (kmaj > 1)
         fprintf(fp, "void ATL_%cu%dam2rm_a%c_%dx%d%s\n",
                 pre, UID, alp, kmaj, u, cst[g]);
      else
         fprintf(fp, "void ATL_%cu%dam2rm_a%c_%d%s\n",pre, UID, alp, u, cst[g]);
      if (pre == 'z' || pre == 'c')
         fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*,const TYPE*);\n");
      else
         fprintf(fp,
         "   (ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*);\n");
   }
}

void SpewForthCpProto(char pre, FILE *fp, char alp, int u, int kmaj)
{
   const int G = (pre == 'c' || pre == 'z') ? 2 : 1;
   const char *cst[2] = {"", "C"};
   int g;

   for (g=0; g < G; g++)
   {
      if (kmaj > 1)
         fprintf(fp, "void ATL_%cu%dcm2am_a%c_%dx%d%s\n",
                 pre, UID, alp, kmaj, u, cst[g]);
      else
         fprintf(fp, "void ATL_%cu%dcm2am_a%c_%d%s\n",pre, UID, alp, u, cst[g]);
      if (pre == 'z' || pre == 'c')
         fprintf(fp,
    "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*,TYPE*);\n");
      else
         fprintf(fp,
         "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*);\n");
      if (kmaj > 1)
         fprintf(fp, "void ATL_%cu%drm2am_a%c_%dx%d%s\n",
                 pre, UID, alp, kmaj, u, cst[g]);
      else
         fprintf(fp, "void ATL_%cu%drm2am_a%c_%d%s\n",pre, UID, alp, u, cst[g]);
      if (pre == 'z' || pre == 'c')
         fprintf(fp,
     "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*,TYPE*);\n");
      else
         fprintf(fp,
         "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*);\n");
   }
}

void SpewForthCpConjDecl(char pre, int REVERSE, ATL_mmnode_t *mmb, FILE *fp,
                         char *arr, char *rt, char alp, int u)
{
   ATL_mmnode_t *mp;
   int j;

   fprintf(fp, "static const %s_t %s_a%c[%d] =\n{\n",
           REVERSE?"am2cm":"cm2am", arr, alp,  ATL_CountNumberOfMMNodes(mmb));
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      const int kmaj = mp->kmaj;
      if (kmaj > 1)
         fprintf(fp, "   ATL_%cu%d%s_a%c_%dx%dC", pre, UID, rt, alp, kmaj,
                 u?mp->mu:mp->nu);
      else
         fprintf(fp, "   ATL_%cu%d%s_a%c_%dC", pre, UID, rt, alp,
                 u?mp->mu:mp->nu);
      if (mp->next)
         fprintf(fp, ",");
      else
         fprintf(fp, " ");
      fprintf(fp, "  /* index %d */\n", j);
   }
   fprintf(fp, "};\n\n");
}

void SpewForthCpDecl(char pre, int REVERSE, ATL_mmnode_t *mmb, FILE *fp,
                     char *arr, char *rt, char alp, int u)
{
   ATL_mmnode_t *mp;
   int j;

   fprintf(fp, "static const %s_t %s_a%c[%d] =\n{\n",
           REVERSE?"am2cm":"cm2am", arr, alp, ATL_CountNumberOfMMNodes(mmb));
   for (j=0,mp=mmb; mp; j++,mp = mp->next)
   {
      const int kmaj = mp->kmaj;
      if (kmaj > 1)
         fprintf(fp, "   ATL_%cu%d%s_a%c_%dx%d", pre, UID, rt, alp, kmaj,
                 u?mp->mu:mp->nu);
      else
         fprintf(fp, "   ATL_%cu%d%s_a%c_%d", pre, UID, rt, alp,
                 u?mp->mu:mp->nu);
      if (mp->next)
         fprintf(fp, ",");
      else
         fprintf(fp, " ");
      fprintf(fp, "  /* index %d */\n", j);
   }
   fprintf(fp, "};\n\n");
}


void GenAMAJ2CMAJ(char pre, ATL_mmnode_t *mmb, char *outd, char *suff)
/*
 * 3. atlas_<pre>amm_am2cm_a[1,X,n]:
 *    defines: ATL_AMM_NCASES
 *    prototypes all am2rm & am2cm routines
 *    1 indexible array giving which to use for each block factor
 */
{
   char ac[3] = {'1', 'n', 'X'};
   int ia, j;
   char *fnam, *sp, *np;
   ATL_mmnode_t *mp;

   if (!suff)
      suff = "";
   ia = strlen(outd) + strlen(suff) + 24+UIL;
   fnam = malloc(ia*sizeof(char));
   assert(fnam);
   sprintf(fnam, "%s/atlas_%cu%damm%s_am2cm_a1.h", outd, pre, UID, suff);
   np = fnam+ia-23+12;
   assert(*np == 'a' && np[1] == 'm');
   sp = fnam+ia-4;
   assert(*sp == '1');

   for (ia=0; ia < 3; ia++)
   {
      char *rt[2] = {"am2cm", "am2rm"};
      FILE *fp;
      int kmaj = mmb->kmaj;

      if (kmaj == 1)
         kmaj = 0;
      *sp = ac[ia];
      fp = fopen(fnam, "w");
      assert(fp);
      sp[1] = '\0';
      PrintBegBlock(pre, mmb, np, fp);
      sp[1] = '.';
      fprintf(fp, "/*\n * mat2blk prototypes\n */\n");
      SpewForthRevCpProto(pre, fp, ac[ia], mmb->mu, kmaj);
      if (mmb->nu != mmb->mu || kmaj > 1)
         SpewForthRevCpProto(pre, fp, ac[ia], mmb->nu, kmaj);
      for (mp=mmb->next; mp; mp = mp->next)
      {
         ATL_mmnode_t *p;
         int mu = mp->mu, nu = mp->nu, kmaj=mp->kmaj;
         for (p=mmb; p != mp; p = p->next)
            if (mu == p->mu && p->kmaj == kmaj)
               break;
         if (p == mp) /* haven't seen before */
            SpewForthRevCpProto(pre, fp, ac[ia], mu, kmaj);
         for (p=mmb; p != mp; p = p->next)
            if ((p->nu == nu && p->kmaj == kmaj) ||
                (p->kmaj == kmaj && p->mu == nu))
               break;
         if (p == mp) /* haven't seen before */
            SpewForthRevCpProto(pre, fp, ac[ia], nu, kmaj);
      }
      fprintf(fp, "\n");
      SpewForthCpDecl(pre,1,mmb, fp, "ATL_UAMM_BLK2A", "am2cm", ac[ia], 1);
      SpewForthCpDecl(pre,1,mmb, fp, "ATL_UAMM_BLK2AT", "am2rm", ac[ia], 1);
      SpewForthCpDecl(pre,1,mmb, fp, "ATL_UAMM_BLK2B", "am2cm", ac[ia], 0);
      SpewForthCpDecl(pre,1,mmb, fp, "ATL_UAMM_BLK2BT", "am2rm", ac[ia], 0);
      if (pre == 'c' || pre == 'z')
      {
         SpewForthCpConjDecl(pre, 1, mmb, fp, "ATL_UAMM_BLKC2A",
                             "am2cm", ac[ia], 1);
         SpewForthCpConjDecl(pre, 1, mmb, fp, "ATL_UAMM_BLKH2A",
                             "am2rm", ac[ia], 1);
         SpewForthCpConjDecl(pre, 1, mmb, fp, "ATL_UAMM_BLKC2B",
                             "am2cm", ac[ia], 0);
         SpewForthCpConjDecl(pre, 1, mmb, fp, "ATL_UAMM_BLKH2B",
                             "am2rm", ac[ia], 0);
      }
      fprintf(fp, "\n#endif  /* end include file guard */\n");
      fclose(fp);
   }
   free(fnam);
}
void GenCMAJ2AMAJ(char pre, ATL_mmnode_t *mmb, char *outd, char *suff)
/*
 * 3. atlas_<pre>amm_cm2am_a[1,X,n]:
 *    defines: ATL_AMM_NCASES
 *    prototypes all rm2am & cm2am routines
 *    1 indexible array giving which to use for each block factor
 */
{
   char ac[3] = {'1', 'n', 'X'};
   int ia, j;
   char *fnam, *sp, *np;
   ATL_mmnode_t *mp;

GenAMAJ2CMAJ(pre, mmb, outd, suff);
   if (!suff)
      suff = "";
   ia = strlen(outd) + strlen(suff) + 24+UIL;
   fnam = malloc(ia*sizeof(char));
   assert(fnam);
   sprintf(fnam, "%s/atlas_%cu%damm%s_cm2am_a1.h", outd, pre, UID, suff);
   np = fnam+ia-23+12;
   assert(*np == 'c' && np[1] == 'm');
   sp = fnam+ia-4;
   assert(*sp == '1');

   for (ia=0; ia < 3; ia++)
   {
      char *rt[2] = {"cm2am", "rm2am"};
      FILE *fp;
      int kmaj = mmb->kmaj;

      if (kmaj == 1)
         kmaj = 0;
      *sp = ac[ia];
      fp = fopen(fnam, "w");
      assert(fp);
      sp[1] = '\0';
      PrintBegBlock(pre, mmb, np, fp);
      sp[1] = '.';
      fprintf(fp, "/*\n * mat2blk prototypes\n */\n");
      SpewForthCpProto(pre, fp, ac[ia], mmb->mu, kmaj);
      if (mmb->nu != mmb->mu || kmaj > 1)
         SpewForthCpProto(pre, fp, ac[ia], mmb->nu, kmaj);
      for (mp=mmb->next; mp; mp = mp->next)
      {
         ATL_mmnode_t *p;
         int mu = mp->mu, nu = mp->nu, kmaj=mp->kmaj;
         for (p=mmb; p != mp; p = p->next)
            if (mu == p->mu && p->kmaj == kmaj)
               break;
         if (p == mp) /* haven't seen before */
            SpewForthCpProto(pre, fp, ac[ia], mu, kmaj);
         for (p=mmb; p != mp; p = p->next)
            if ((p->nu == nu && p->kmaj == kmaj) ||
                (p->kmaj == kmaj && p->mu == nu))
               break;
         if (p == mp) /* haven't seen before */
            SpewForthCpProto(pre, fp, ac[ia], nu, kmaj);
      }
      fprintf(fp, "\n");
      SpewForthCpDecl(pre,0,mmb, fp, "ATL_UAMM_A2BLK", "cm2am", ac[ia], 1);
      SpewForthCpDecl(pre,0,mmb, fp, "ATL_UAMM_AT2BLK", "rm2am", ac[ia], 1);
      SpewForthCpDecl(pre,0,mmb, fp, "ATL_UAMM_B2BLK", "cm2am", ac[ia], 0);
      SpewForthCpDecl(pre,0,mmb, fp, "ATL_UAMM_BT2BLK", "rm2am", ac[ia], 0);
      if (pre == 'c' || pre == 'z')
      {
         SpewForthCpConjDecl(pre, 0, mmb, fp, "ATL_UAMM_AC2BLK",
                             "cm2am", ac[ia], 1);
         SpewForthCpConjDecl(pre, 0, mmb, fp, "ATL_UAMM_AH2BLK",
                             "rm2am", ac[ia], 1);
         SpewForthCpConjDecl(pre, 0, mmb, fp, "ATL_UAMM_BC2BLK",
                             "cm2am", ac[ia], 0);
         SpewForthCpConjDecl(pre, 0, mmb, fp, "ATL_UAMM_BH2BLK",
                             "rm2am", ac[ia], 0);
      }
      fprintf(fp, "\n#endif  /* end include file guard */\n");
      fclose(fp);
   }
   free(fnam);
}

int KernelIsExactSame(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if kernels are the same including KB, 0 otherwise
 */
{
/*
 * Kernels aren't the same if one is being compiled with specific KB,
 * and the other has runtime
 */
   if (FLAG_IS_SET(p0->flag, MMF_KRUNTIME) !=
       FLAG_IS_SET(p1->flag, MMF_KRUNTIME))
      return(0);
   if (!FLAG_IS_SET(p0->flag, MMF_KRUNTIME) && (p0->kbB != p1->kbB))
      return(0);
/*
 * Kernels aren't same if one is K-veclen-major and one is not
 */
   if (p0->kmaj != p1->kmaj)
      return(0);
/*
 * Kernels aren't same if the -DATL_MOVE bits don't match
 */
   if (ATL_MMF_MVGET(p0->flag) != ATL_MMF_MVGET(p1->flag))
      return(0);
/*
 * Two generated kernels are the same if mu,nu,kmaj,ku,VLEN,flag are the same.
 * NOTE: if we make generator handle muladd, etc, MUST UPDATE HERE!!!
 */
   if (p0->ID == 0 && p1->ID == 0)
      return(p0->mu == p1->mu && p0->nu == p1->nu && p0->kmaj == p1->kmaj &&
             p0->ku == p1->ku && p0->vlen == p1->vlen && p0->flag == p1->flag);
/*
 * If both are user kernels, then they may be repeats.  For user kernels,
 * they are the same if both ID and flag match, else they are not.
 */
   else if (p0->ID > 0 && p1->ID > 0)
      return(p0->ID == p1->ID && p0->flag == p1->flag);
   return(0);  /* Can't be the same if above criteria fails */
}

/*
 * RETURNS: flags that necessitate recompilation, not including KRUNTIME,
 * which is encoded in kb
 */
int GetCompTimeFlags(ATL_mmnode_t *mp)
{
   int iflg;
   iflg = ATL_MMF_MVGET(mp->flag);  /* MVbits change kern at comp time */
   iflg |=  (((mp->flag) & 1)<<3);  /* LDTOP/BOT could be compile-time dec */
   if (mp->kmaj > 1)
      iflg |= 1<<4;
   return(iflg);
}
int ExactKernelInList(ATL_mmnode_t *mmb, ATL_mmnode_t *p)
/*
 * RETURNS: 1 if p is duplicated in mmb, else 0
 */
{
   ATL_mmnode_t *mp;
   if (!p || !mmb)
      return(0);
   for (mp=mmb; mp; mp = mp->next)
      if (KernelIsExactSame(mp, p))
         return(1);
    return(0);
}

void SpewForthKernProto(FILE *fp, char pre, ATL_mmnode_t *p, char bc)
{
   fprintf(fp, "void ATL_%cu%dAMMM_%d_%d_%x_%dx%dx%d_b%c\n", pre, UID, p->ID,
           FLAG_IS_SET(p->flag, MMF_KRUNTIME)?0:p->kbB, GetCompTimeFlags(p),
           p->mu, p->nu, p->ku,bc);
   fprintf(fp,
      "   (ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,TYPE*,\n");
   fprintf(fp,
      "    const TYPE*,const TYPE*,const TYPE*);\n");
}

void SpewForthKernProtos(FILE *fp, char pre, ATL_mmnode_t *mmb, int nbet)
{
   ATL_mmnode_t *mp;
   for (mp=mmb; mp; mp = mp->next)
   {
      if (!ExactKernelInList(mp->next, mp))
      {
         char bc[3] = {'0', '1', 'n'};  /* 0 must come first */
         int ib;
         for (ib=0; ib < nbet; ib++)
            SpewForthKernProto(fp, pre, mp, bc[ib]);
      }
   }
}

void SpewForthKernArray(FILE *fp, char pre, ATL_mmnode_t *mmb,
                        char *vnam, char cbet)
{
   ATL_mmnode_t *mp;
   int n;

   for (n=0,mp=mmb; mp; n++, mp = mp->next);
   fprintf(fp, "static const ammkern_t ATL_UAMM_%s[%d] =\n", vnam, n);
   fprintf(fp, "{\n");
   for (mp=mmb; mp; mp = mp->next)
   {
      fprintf(fp, "   ATL_%cu%dAMMM_%d_%d_%x_%dx%dx%d_b%c", pre, UID, mp->ID,
              FLAG_IS_SET(mp->flag, MMF_KRUNTIME)?0:mp->kbB,
              GetCompTimeFlags(mp), mp->mu, mp->nu, mp->ku, cbet);
      if (mp->next)
         fprintf(fp, ",\n");
   }
   fprintf(fp, "\n};\n\n");
}

/*
 * RETURNS: possibly updated list of all unique mu/nu comboes
 */
typedef struct mnur mnur_t;
struct mnur {int mu; int nu; int kmaj; mnur_t *next;};
mnur_t *GetUniqueMNUnrolls(ATL_mmnode_t *mmb, mnur_t *urb)
{
   ATL_mmnode_t *mp;
   mnur_t *up;
/*
 * For each node in mmb, add to urb if mu/nu combo not already there
 * kmaj only affects A/B copy, and this is for C put, so ignore kmaj
 */
   for (mp = mmb; mp; mp = mp->next)
   {
      for (up=urb; up; up = up->next)
         if (mp->mu == up->mu && mp->nu == up->nu)
            break;
      if (!up)
      {
         up = malloc(sizeof(mnur_t));
         assert(up);
         up->mu = mp->mu;
         up->nu = mp->nu;
         up->kmaj = 0;
         up->next = urb;
         urb = up;
      }
   }
   return(urb);
}

/*
 * RETURNS: list of just unique MUs from mnb not already in mub
 */
mnur_t *GetUniqueMUnrolls(ATL_mmnode_t *mnb, mnur_t *mub)
{
   mnur_t *mup;
   ATL_mmnode_t *mnp;
   if (!mnb)
      return(mub);
   for (mnp=mnb; mnp; mnp = mnp->next)
   {
      const int kmaj = mnp->kmaj;
      for (mup=mub; mup; mup = mup->next)
         if (mup->mu == mnp->mu && mup->kmaj == kmaj)
            break;
      if (!mup)  /* a new mu */
      {
         mup = malloc(sizeof(mnur_t));
         assert(mup);
         mup->nu = mup->mu = mnp->mu;
         mup->next = mub;
         mup->kmaj = kmaj;
         mub = mup;
      }
   }
   return(mub);
}
/*
 * RETURNS: list of just unique NUs from mnb not already in mub
 */
mnur_t *GetUniqueNUnrolls(ATL_mmnode_t *mnb, mnur_t *mub)
{
   mnur_t *mup;
   ATL_mmnode_t *mnp;
   if (!mnb)
      return(mub);
   for (mnp=mnb; mnp; mnp = mnp->next)
   {
      const int kmaj = mnp->kmaj;
      for (mup=mub; mup; mup = mup->next)
         if (mup->mu == mnp->nu && kmaj == mup->kmaj)
            break;
      if (!mup)  /* a new nu */
      {
         mup = malloc(sizeof(mnur_t));
         assert(mup);
         mup->nu = mup->mu = mnp->nu;
         mup->next = mub;
         mup->kmaj = kmaj;
         mub = mup;
      }
   }
   return(mub);
}

void KillUnrollList(mnur_t *b)
{
   mnur_t *p;
   while (b)
   {
      p = b->next;
      free(b);
      b = p;
   }
}
void PrintSwapProto(FILE *fp, char pre, int mu, int nu)
{
   fprintf(fp, "void Mjoin(PATL,ammswp%dx%d)", mu, nu);
   if (pre == 'd' || pre == 's')
      fprintf(fp, "(ATL_CINT nnu, TYPE *A, ATL_CSZT lda, TYPE *b);\n");
   else
      fprintf(fp, "(ATL_CINT nnu, TYPE *A, ATL_CSZT lda, TYPE *r, TYPE *i);\n");
}


void zGenAmmSwp(char pre, FILE *fp, int mu, int nu)
{
}

void GenAmmSwp(char pre, FILE *fp, int mu, int nu)
{
   const int munu=mu*nu;
   int i, j, ib;

   if (pre == 'z' || pre == 'c')
   {
      zGenAmmSwp(pre, fp, mu, nu);
      return;
   }

   fprintf(fp, "#include \"atlas_misc.h\"\n");
   fprintf(fp, "void Mjoin(PATL,ammswp%dx%d)\n", mu, nu);
   fprintf(fp, "(\n");
   fprintf(fp, "   ATL_CINT nnu,   /* CEIL(rowlen / nu) */\n");
   fprintf(fp, "   TYPE *A,        /* col-maj matrix to swap wt b */\n");
   fprintf(fp, "   ATL_CSZT lda1,  /* stride between row elts in A */\n");
   fprintf(fp,
           "   TYPE *b         /* %dx%d C-format row ptr to be swapped */\n",
           mu, nu);
   fprintf(fp, ")\n{\n");

   fprintf(fp, "   register unsigned int j;\n");
   if (nu > 1)
   {
      fprintf(fp, "   const size_t lda2=lda1+lda1");
      for (i=3; i <= nu; i++)
         fprintf(fp, ", lda%d=lda1+lda%d", i, i-1);
      fprintf(fp, ";\n");
   }

   fprintf(fp, "   for (j=nnu; j; j--)\n   {\n");

   fprintf(fp, "      register TYPE a0");
   for (i=1; i < nu; i++)
      fprintf(fp, ", a%d", i);
   fprintf(fp, ";\n");

   fprintf(fp, "      a0 = *A;\n");
   for (i=1; i < nu; i++)
      fprintf(fp, "      a%d = A[lda%d];\n", i, i);
   fprintf(fp, "      *A = *b;\n");
   for (i=1; i < nu; i++)
      fprintf(fp, "      A[lda%d] = b[%d];\n", i, i*mu);
   fprintf(fp, "      *b = a0;\n");
   for (i=1; i < nu; i++)
      fprintf(fp, "      b[%d] = a%d;\n", i*mu, i);
   fprintf(fp, "      b += %d;\n", mu*nu);
   fprintf(fp, "      A += lda%d;\n", nu);

   fprintf(fp, "   }\n");

   fprintf(fp, "}\n");
}

void GenAmmSwapFiles(char pre, ATL_mmnode_t *mmb, char *outd)
{
   mnur_t *ub, *up;
   char *fnam;
   int ia;

   if (pre == 's')
      GenAmmSwapFiles('c', mmb, outd);
   else if (pre == 'd')
      GenAmmSwapFiles('z', mmb, outd);
   ia = strlen(outd) + 24+UIL;
   fnam = malloc(ia);
   assert(fnam);
   ub = GetUniqueMNUnrolls(mmb, NULL);
   for (up=ub; up; up = up->next)
   {
      FILE *fp;
      assert(up->mu < 100 && up->nu < 100);
      sprintf(fnam, "%s/ATL_%cammswp%dx%d.c", outd, pre, up->mu, up->nu);
      fp = fopen(fnam, "w");
      assert(fp);
      GenAmmSwp(pre, fp, up->mu, up->nu);
      fclose(fp);
   }
   KillUnrollList(ub);
   free(fnam);
}
void GenAmmSwapH(char pre, ATL_mmnode_t *mmb, char *outd)
{
   mnur_t *ub, *up;
   ATL_mmnode_t *mp;
   FILE *fp;
   int i;

   fp = StandHStart(pre, mmb, outd, "swp");

   fprintf(fp, "\n");
   ub = GetUniqueMNUnrolls(mmb, NULL);
   for (up=ub; up; up = up->next)
      PrintSwapProto(fp, pre, up->mu, up->nu);
   KillUnrollList(ub);
   fprintf(fp, "\n");

   fprintf(fp, "static const ammswp_t ATL_UAMM_SWP[ATL_UAMM_NCASES] =\n{\n");
   for (i=0, mp=mmb; mp; mp = mp->next, i++)
      if (mp->next)
         fprintf(fp, "   ATL_%cammswp%dx%d,  /* index %d */\n",
                 pre, mp->mu, mp->nu, i);
     else
         fprintf(fp, "   ATL_%cammswp%dx%d   /* index %d */\n",
                 pre, mp->mu, mp->nu, i);

   fprintf(fp, "};\n\n");
   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

int KernelIsSame(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if kernels are the same except for blocking, 0 otherwise
 */
{
/*
 * Kernels aren't the same if one is being compiled with specific KB,
 * and the other has runtime
 */
   if (FLAG_IS_SET(p0->flag, MMF_KRUNTIME) !=
       FLAG_IS_SET(p1->flag, MMF_KRUNTIME))
      return(0);
/*
 * Two generated kernels are the same if mu,nu,ku,VLEN,flag are the same.
 * NOTE: if we make generator handle muladd, etc, MUST UPDATE HERE!!!
 */
   if (p0->ID == 0 && p1->ID == 0)
      return(p0->mu == p1->mu && p0->nu == p1->nu && p0->ku == p1->ku &&
             p0->vlen == p1->vlen && p0->flag == p1->flag);
/*
 * If both are user kernels, then they may be repeats.  For user kernels,
 * they are the same if both ID and flag match, else they are not.
 */
   else if (p0->ID > 0 && p1->ID > 0)
      return(p0->ID == p1->ID && p0->flag == p1->flag);
   return(0);  /* Can't be the same if above criteria fails */
}


void GenRankKH
(
   char pre,
   ATL_mmnode_t *sqb,  /* baseptr  of square-case AMMM kernels */
   ATL_mmnode_t *rkb,  /* rank-K kernels, one for each supported K */
   char *outd
)
{
   FILE *fp;
   mnur_t *putb, *cpyb, *up;
   ATL_mmnode_t *mp;
   int m, n, k;
   int ia, ib;
   char PRE = pre;
   char ac[3] =  {'1', 'n', 'X'};
   char bc[4] = {'1', 'n', '0', 'X'};
   if (pre == 'c')
      pre = 's';
   else  if (pre == 'z')
      pre = 'd';

   assert(sqb && rkb);
   fp = StandHStart(PRE, rkb, outd, "rankK");
   fprintf(fp, "\n");

   for (mp=sqb; mp->next; mp = mp->next);
   fprintf(fp, "#define ATL_MAXM_RKK %d\n", mp->mbB);
   fprintf(fp, "#define ATL_MAXN_RKK %d\n", mp->nbB);
   fprintf(fp, "#define ATL_MAXK_RKK %d\n\n", mp->kbB);
/*
 * Prototype needed copy routines
 */
   putb = GetUniqueMNUnrolls(rkb, NULL);
   fprintf(fp, "/*\n * cblk2mat put function prototypes\n */\n");
   for (up=putb; up; up = up->next)
   {
      for (ib=0; ib < 4; ib++)
      {
         fprintf(fp, "void ATL_%cablk2cmat_%dx%d_a1_b%c\n",
                 PRE, up->mu, up->nu, bc[ib]);
         if (PRE == 'c' || PRE == 'z')
            fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const TYPE*,const SCALAR,TYPE *,ATL_CSZT);\n");
         else
            fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const SCALAR,TYPE *,ATL_CSZT);\n");
      }
   }
   KillUnrollList(putb);
   fprintf(fp,
      "/*\n * Column-major to access-major copy function prototypes\n */\n");
   cpyb = GetUniqueMUnrolls(rkb, NULL);
   cpyb = GetUniqueNUnrolls(rkb, cpyb);
   for (up=cpyb; up; up = up->next)
   {
      for (ia=0; ia < 3; ia++)
         SpewForthCpProto(PRE, fp, ac[ia], up->mu, up->kmaj);
   }
   KillUnrollList(cpyb);
/*
 * Prototype the rank-K functions
 */
   fprintf(fp, "/*\n * rank-K AMMM kernel prototypes\n */\n");
   for (mp=rkb; mp; mp = mp->next)
   {
/*
 *    For runtime kernels, only prototype 1st time they are seen in list
 */
      if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
      {
         ATL_mmnode_t *p;
         if (mp->ID > 0)
         {
            for (p=rkb; p != mp; p = p->next)
               if (p->ID == mp->ID)
                  break;
         }
         else
         {
            for (p=rkb; p != mp; p = p->next)
               if (KernelIsSame(mp, p))
                  break;
         }
         if (mp == p)
         {
            SpewForthKernProto(fp, pre, mp, '0');
            SpewForthKernProto(fp, pre, mp, '1');
            SpewForthKernProto(fp, pre, mp, 'n');
         }
      }
/*
 *    compile-time-K kernels get prototyped for each invocation
 */
      else
      {
         SpewForthKernProto(fp, pre, mp, '0');
         SpewForthKernProto(fp, pre, mp, '1');
         SpewForthKernProto(fp, pre, mp, 'n');
      }
   }
/*
 * Now, crank out funcptr arrays
 */
   for (ib=0; ib < 4; ib++)
   {
      fprintf(fp,
         "\nstatic const ablk2cmat_t ATL_RKK_BLK2C_a1_b%c[%d] =\n{\n",
              bc[ib], ATL_CountNumberOfMMNodes(rkb));
      SpewForthC2BDecl(PRE, rkb, fp, "ablk2cmat", '1', bc[ib]);
   }
   for (ia=0; ia < 3; ia++)
   {
      fprintf(fp, "\n");
      SpewForthCpDecl(PRE, 0, rkb, fp, "ATL_RKK_A2BLK", "cm2am", ac[ia], 1);
      SpewForthCpDecl(PRE, 0, rkb, fp, "ATL_RKK_AT2BLK", "rm2am", ac[ia], 1);
      SpewForthCpDecl(PRE, 0, rkb, fp, "ATL_RKK_B2BLK", "cm2am", ac[ia], 0);
      SpewForthCpDecl(PRE, 0, rkb, fp, "ATL_RKK_BT2BLK", "rm2am", ac[ia], 0);
      if (PRE == 'z' || PRE == 'c')
      {
         SpewForthCpConjDecl(PRE,0,rkb, fp, "ATL_RKK_AC2BLK", "cm2am",ac[ia],1);
         SpewForthCpConjDecl(PRE,0,rkb, fp, "ATL_RKK_AH2BLK", "rm2am",ac[ia],1);
         SpewForthCpConjDecl(PRE,0,rkb, fp, "ATL_RKK_BC2BLK", "cm2am",ac[ia],0);
         SpewForthCpConjDecl(PRE,0,rkb, fp, "ATL_RKK_BH2BLK", "rm2am",ac[ia],0);
      }
   }
   SpewForthKernArray(fp, pre, rkb, "KERN_RKK", '0');
   SpewForthKernArray(fp, pre, rkb, "KERN_RKK_b1", '1');
   SpewForthKernArray(fp, pre, rkb, "KERN_RKK_bn", 'n');
   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenSquareKH(char pre, ATL_mmnode_t *mmb, char *outd)
/*
 * 5. atlas_<pre>amm_kerns.h
 *    defines: ATL_AMM_NCASES
 *    prototypes all kernels, including K-cleanup
 *    1 indexible array gives kernel to use for each case as func ptr
 *    1 indexible array gives K-clean kernel
 */
{
   FILE *fp;
/*
 * Dump out standard header start and kernel prototypes
 */
   fp = StandHStart(pre, mmb, outd, "sqkern");
   fprintf(fp, "\n");
   SpewForthKernProtos(fp, pre, mmb, 3);
   fprintf(fp, "\n");
/*
 * Dump out kernel ptr arrays
 */
   SpewForthKernArray(fp, pre, mmb, "SQKERN_b0", '0');
   SpewForthKernArray(fp, pre, mmb, "SQKERN_b1", '1');
   SpewForthKernArray(fp, pre, mmb, "SQKERN_bn", 'n');

   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenKernH(char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *ukb,
              ATL_mmnode_t *kcb, char *outd)
/*
 * 5. atlas_<pre>amm_kerns.h
 *    defines: ATL_AMM_NCASES
 *    prototypes all kernels, including K-cleanup
 *    1 indexible array gives kernel to use for each case as func ptr
 *    1 indexible array gives K-clean kernel
 */
{
   FILE *fp;
/*
 * Dump out standard header start and kernel prototypes
 */
   fp = StandHStart(pre, mmb, outd, "kern");
   fprintf(fp, "\n");
   SpewForthKernProtos(fp, pre, mmb, 3);
   if (ukb)
      SpewForthKernProtos(fp, pre, ukb, 3);
   fprintf(fp, "\n");
/*
 * Dump out kernel ptr arrays
 */
   SpewForthKernArray(fp, pre, mmb, "KERN_b0", '0');
   SpewForthKernArray(fp, pre, mmb, "KERN_b1", '1');
   SpewForthKernArray(fp, pre, mmb, "KERN_bn", 'n');
   if (kcb)
   {
      SpewForthKernArray(fp, pre, kcb, "KERN_K1", '0');
      SpewForthKernArray(fp, pre, kcb, "KERN_K1_b1", '1');
      SpewForthKernArray(fp, pre, kcb, "KERN_K1_bn", 'n');
   }

   fprintf(fp, "\n#endif  /* end include file guard */\n");
   fclose(fp);
}

void GenCmplxHeaders(char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *rkb, char *outd)
{
   GenCMAJ2AMAJ(pre, mmb, outd, NULL);
   GenC2BLK(pre, mmb, outd, NULL);
   if (rkb)
   {
      GenCMAJ2AMAJ(pre, rkb, outd, "rkk");
      GenC2BLK(pre, rkb, outd, "_rkk");
      GenRankKH(pre, mmb, rkb, outd);
   }
}

void GenHeaderFiles(char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *ukb,
                    ATL_mmnode_t *kcb, ATL_mmnode_t *rkb, ATL_mmnode_t *urb,
                    ATL_mmnode_t *sqb, ATL_mmnode_t *usb, char *outd)
/*
 * Header files required to build full gemm (no timing):
 *X1. atlas_<pre>amm_blk.h :
 *X   defines: ATL_AMM_NCASES, ATL_AMM_MAX[M,N,K]B
 *X   3 arrays indexed by case give blocking
 *X2  atlas_<pre>amm_flag.h
 *X   defines: ATL_AMM_NCASES
 *X   1 indexible array giving KRUNTIME for now
 *X3. atlas_<pre>amm_cm2am_a[1,X,n]:
 *X   defines: ATL_AMM_NCASES
 *X   prototypes all rm2am & cm2am routines
 *X   1 indexible array giving which to use for each block factor
 *X4. atlas_<pre>amm_ablk2cmat.h
 *X   defines: ATL_AMM_NCASES
 *X   prototypes all ablk2cmat routines
 *X   1 indexible array for each alpha,beta combination
 *X   -> 3*4 = 12 indexible arrays total
 *X5. atlas_<pre>amm_kerns.h
 *X   defines: ATL_AMM_NCASES
 *X   prototypes all kernels, including K-cleanup
 *X   1 indexible array gives kernel to use for each case as func ptr
 *X   1 indexible array gives K-clean kernel
 *X6. atlas_<pre>amm_cmat2ablk.h (I don't need, Rakib does)
 *X   defines: ATL_AMM_NCASES
 *X   prototypes all cmat2ablk routines
 *X   1 indexible array for each alpha,beta combination
 *X   -> 3*4 = 12 indexible arrays total
 */
{
   if (rkb)
      GenAmmSum(pre, mmb, rkb, outd);
   GenAmmSwapH(pre, mmb, outd);
   GenBlockingFile(pre, mmb, outd, "blk");
   GenPerfFile(pre, mmb, outd, "perf");
   GenFlagH(pre, mmb, outd, "flag");
   GenCMAJ2AMAJ(pre, mmb, outd, NULL);
   GenC2BLK(pre, mmb, outd, NULL);
   GenKernH(pre, mmb, ukb, kcb, outd);
   if (rkb)
   {
      GenBlockingFile(pre, rkb, outd, "rkkblk");
      GenFlagH(pre, rkb, outd, "rkkflag");
      GenCMAJ2AMAJ(pre, rkb, outd, "rkk");
      GenC2BLK(pre, rkb, outd, "_rkk");
      GenRankKH(pre, mmb, rkb, outd);
   }
   if (sqb)
   {
      GenBlockingFile(pre, rkb, outd, "sqblk");
      GenFlagH(pre, rkb, outd, "sqflag");
      GenCMAJ2AMAJ(pre, rkb, outd, "sq");
      GenC2BLK(pre, rkb, outd, "_sq");
      GenSquareKH(pre, sqb, outd);
   }
}

/*
 * Splits rkb into two lists: (1) Routines with runtime K (RUNB),
 * (2) Routines with fixed-KB (returned)
 * NOTE: original list rkb is unchanged
 */
ATL_mmnode_t *SplitRankK(ATL_mmnode_t *rkb, ATL_mmnode_t **RUNB)
{
   ATL_mmnode_t *runb=NULL, *fixb=NULL, *p, *np;
   for (p=rkb; p; p = p->next)
   {
      np = CloneMMNode(p);
      if (FLAG_IS_SET(np->flag, MMF_KRUNTIME))
      {
         np->next = runb;
         runb = np;
      }
      else
      {
         np->next = fixb;
         fixb = np;
      }
   }
   *RUNB = runb;
   return(fixb);
}

char *GetKernComp(ATL_mmnode_t *mmp, char *dcomp, char *dflags, char **flgs)
{
   char *comp = dcomp;
   if (mmp->comp)
   {
      comp = (mmp->comp[0] == 'g' && mmp->comp[1] == 'c' &&
              mmp->comp[2] == 'c' &&
             (mmp->comp[3] == '\0' || mmp->comp[3] == ' '))
             ? "$(GOODGCC)" : mmp->comp;
      *flgs = mmp->cflags;
   }
   else
      *flgs = dflags;
   return(comp);
}
void PrintKernComp
(
   FILE *fp,            /* file to print to */
   char pre,
   ATL_mmnode_t *mmp,   /* kernel compile rule is for */
   int UID,             /* user-ID for user-determined kerns */
   char *comp,          /* compiler to use */
   char *cflags,        /* compiler flags to use */
   char *styp,          /* string defining type (eg. "-DSREAL") */
   char cbet,           /* character with beta name ('1', '0', 'n') */
   char *sbet           /* string wt full beta name ("1", "N1", "0") */
)
{
   const int kb = FLAG_IS_SET(mmp->flag, MMF_KRUNTIME)?0:mmp->kbB;
   const int flg = GetCompTimeFlags(mmp);
   if (pre == 'z')
      styp = "-DDREAL=1";
   else if (pre == 'c')
      styp = "-DSREAL=1";
   fprintf(fp, "ATL_%cu%dAMMM_%d_%d_%x_%dx%dx%d_b%c.o : %s\n", pre,
           UID, mmp->ID, kb, flg, mmp->mu, mmp->nu, mmp->ku, cbet, mmp->rout);
   fprintf(fp, "\t%s $(CDEFS2) %s -DBETA%s=1", comp, styp, sbet);
   if (!FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
      fprintf(fp, " -DMB=%d -DNB=%d, -DKB=%d", mmp->mbB, mmp->nbB, mmp->kbB);
   if (FLAG_IS_SET(mmp->flag, MMF_MVA))
      fprintf(fp, " -DATL_MOVEA");
   if (FLAG_IS_SET(mmp->flag, MMF_MVB))
      fprintf(fp, " -DATL_MOVEB");
   if (FLAG_IS_SET(mmp->flag, MMF_MVC))
      fprintf(fp, " -DATL_MOVEC");
         fprintf(fp,
         " -DATL_USERMM=ATL_%cu%dAMMM_%d_%d_%x_%dx%dx%d_b%c -DATL_UAMMID=%d",
                 pre, UID, mmp->ID, kb, flg, mmp->mu, mmp->nu, mmp->ku, cbet,
                 UID);
         fprintf(fp, " %s -o ATL_%cu%dAMMM_%d_%d_%x_%dx%dx%d_b%c.o -c %s\n",
                 cflags, pre, UID, mmp->ID, kb, flg, mmp->mu, mmp->nu,
                 mmp->ku, cbet, mmp->rout);
}
void GenMakefile
(
   char pre,            /* type/precision prefix : s,d,c,z */
   ATL_mmnode_t *mmb,   /* main kernels for GEMM */
   ATL_mmnode_t *ukb,   /* unique kernels for doing K cleanup of kerns in mmb */
   ATL_mmnode_t *rkb,   /* list of kernels to doing rank-K update */
   ATL_mmnode_t *urb,   /* rank-K update kerns not existing in other lists */
   ATL_mmnode_t *usb,   /* square kernels not existing in other lists */
   char *outd
)
{
   ATL_mmnode_t *mmp, *p, *fixb, *runb;
   mnur_t *mnurb=NULL, *allub, *up;
   FILE *fp;
   char *comp, *cflags;
   char *ln;
   int i;
   char pres[2];
   char be[3] = {'1', 'n', '0'};
   char *bes[3] = {"1", "N1", "0"};
   char al[3] = {'1', 'n', 'X'};
   char dcomp[8] = {'$', '(', 'D', 'M', 'C', ')', '\0'};
   char dflags[12] = {'$', '(', 'D', 'M', 'C', 'F', 'L', 'A', 'G', 'S',
                     ')', '\0'};
   char *styps[2] = {"-DDREAL", "-DDCPLX"};
   char *styp = (pre == 'd' || pre == 'z') ? "-DDREAL" : "-DSREAL";

   pres[0] = pre;
   if (pre == 's' || pre == 'c')
   {
      styps[0] = "-DSREAL";
      styps[1] = "-DSCPLX";
      pres[1] = 'c';
   }
   else
      pres[1] = 'z';

   ln = malloc((strlen(outd)+11)*sizeof(char));
   assert(ln);
   sprintf(ln, "%s/%cMake_amm", outd, pre);
   fp = fopen(ln, "w");
   assert(fp);
   free(ln);
   fprintf(fp, "include ../Make.inc\n");
   fprintf(fp, "CDEFS2=$(CDEFS)\n\n");
   if (pre == 'c')
   {
      fprintf(fp, "CMC=$(SMC)\n");
      fprintf(fp, "CKCFLAGS=$(SKCFLAGS)\n");
      fprintf(fp, "CMCFLAGS=$(SMCFLAGS)\n");
   }
   else if (pre == 'z')
   {
      fprintf(fp, "ZMC=$(DMC)\n");
      fprintf(fp, "ZKCFLAGS=$(DKCFLAGS)\n");
      fprintf(fp, "ZMCFLAGS=$(DMCFLAGS)\n");
   }
/*
 * Build list of all unique MU/NU combos for copy routines
 * Square cases built from mmb, so they are all represented in mmb
 */
   mnurb = GetUniqueMNUnrolls(mmb, NULL);
   mnurb = GetUniqueMNUnrolls(urb, mnurb);
   allub = GetUniqueMUnrolls(mmb, NULL);
   allub = GetUniqueMUnrolls(urb, allub);
   allub = GetUniqueNUnrolls(mmb, allub);
   allub = GetUniqueNUnrolls(urb, allub);
/*
 * Spew out all filenames that must be compiled
 */
   fprintf(fp, "objs =");
/*
 * Routines to copy from MU/NU-major to column major output array
 */
   for (up=mnurb; up; up = up->next)
   {
      int j;
      const int mu=up->mu, nu=up->nu;

      for (j=0; j < 3; j++)
      {
         int k;
         char *rtn[2] = {"ablk2cmat", "cmat2ablk"};
         for (k=0; k < 2; k++)
         {
            int h;
            for (h=0; h < 2; h++)
            {
               if (pre != pres[h])
                  continue;
               if (j == 0)
                  fprintf(fp, " \\\n       ATL_%cammswp%dx%d.o", pre, mu, nu);
               fprintf(fp, " \\\n       ATL_%cu%d%s_%dx%d_a%c_b1.o",
                       pre, UID, rtn[k], mu, nu, al[j]);
               fprintf(fp, " ATL_%cu%d%s_%dx%d_a%c_bX.o",
                       pre, UID, rtn[k], mu, nu, al[j]);
               fprintf(fp, " \\\n       ATL_%cu%d%s_%dx%d_a%c_b0.o",
                       pre, UID, rtn[k], mu, nu, al[j]);
               fprintf(fp, " ATL_%cu%d%s_%dx%d_a%c_bn.o",
                       pre, UID, rtn[k], mu, nu, al[j]);
            }
         }
      }
   }
/*
 * Routines to copy back and forth from A and B
 */
   for (up=allub; up; up = up->next)
   {
      int h;
      for (h=0; h < 2; h++)
      {
        if (pre != pres[h]) continue;
         int j;
         const int u = up->mu;
         for (j=0; j < 3; j++)
         {
            if (up->kmaj > 1)
            {
               fprintf(fp,
 " \\\n       ATL_%cu%drm2am_a%c_%dx%d.o ATL_%cu%dcm2am_a%c_%dx%d.o",
              pre, UID, al[j], up->kmaj, u, pre, UID, al[j], up->kmaj, u);
               fprintf(fp,
 " \\\n       ATL_%cu%dam2rm_a%c_%dx%d.o ATL_%cu%dam2cm_a%c_%dx%d.o",
              pre, UID, al[j], up->kmaj, u, pre, UID, al[j], up->kmaj, u);
            }
            else
            {
               fprintf(fp,
                  " \\\n       ATL_%cu%drm2am_a%c_%d.o ATL_%cu%dcm2am_a%c_%d.o",
                       pre, UID, al[j], u, pre, UID, al[j], u);
               fprintf(fp,
                  " \\\n       ATL_%cu%dam2rm_a%c_%d.o ATL_%cu%dam2cm_a%c_%d.o",
                       pre, UID, al[j], u, pre, UID, al[j], u);
            }
            if (pre == 'c' || pre == 'z')
            {
               if (up->kmaj > 1)
               {
                  fprintf(fp,
 " \\\n       ATL_%cu%drm2am_a%c_%dx%dC.o ATL_%cu%dcm2am_a%c_%dx%dC.o",
                  pre, UID, al[j], up->kmaj, u, pre, UID, al[j], up->kmaj, u);
                  fprintf(fp,
 " \\\n       ATL_%cu%dam2rm_a%c_%dx%dC.o ATL_%cu%dam2cm_a%c_%dx%dC.o",
                  pre, UID, al[j], up->kmaj, u, pre, UID, al[j], up->kmaj, u);
               }
               else
               {
                  fprintf(fp,
               " \\\n       ATL_%cu%drm2am_a%c_%dC.o ATL_%cu%dcm2am_a%c_%dC.o",
                          pre, UID, al[j], u, pre, UID, al[j], u);
                  fprintf(fp,
               " \\\n       ATL_%cu%dam2rm_a%c_%dC.o ATL_%cu%dam2cm_a%c_%dC.o",
                          pre, UID, al[j], u, pre, UID, al[j], u);
               }
            }
         }
      }
   }
/*
 * AMM kernel routines
 */
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      int kb = mmp->kbB;
/*
 *    Kernels that take runtime K are only compiled once, so don't repeat them
 *    for every KB.  Only generate a statement if this is the first one.
 */
      if (FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
      {
         const int id = mmp->ID;
         for (p=mmb; p != mmp; p = p->next)
            if (p->ID == id && FLAG_IS_SET(p->flag, MMF_KRUNTIME))
               break;
         if (p != mmp)
            continue;
         kb = 0;
      }
/*
 *    ATL_<pre>UAMMM_<ID>_<kb>_<flg>_<mu>x<nu>x<ku>_b<X>
 */
      for (i=0; i < 3; i++)
         fprintf(fp, " \\\n       ATL_%cu%dAMMM_%d_%d_%x_%dx%dx%d_b%c.o",
                 pre, UID, mmp->ID, kb, GetCompTimeFlags(mmp),
                 mmp->mu, mmp->nu, mmp->ku, be[i]);
   }
/*
 * AMM K-cleanup kernel routines are all unique, so no checking for repeats
 *    ATL_<pre>UAMMM_<ID>_<kb>_<flg>_<mu>x<nu>x<ku>_b0
 */
   for (mmp=ukb; mmp; mmp = mmp->next)
      for (i=0; i < 3; i++)
         fprintf(fp, " \\\n       ATL_%cu%dAMMM_%d_0_%x_%dx%dx%d_b%c.o", pre,
                 UID, mmp->ID, GetCompTimeFlags(mmp), mmp->mu, mmp->nu,
                 mmp->ku, be[i]);
/*
 * library make targets
 */
   fprintf(fp, "\n\nlib : %clib.grd\nall : %clib.grd\n%clib : %clib.grd\n",
           pre, pre, pre, pre);
   fprintf(fp, "%clib.grd : $(objs)\n", pre);
   fprintf(fp, "\t$(ARCHIVER) $(ARFLAGS) $(UAMMlib) $(objs)\n");
   fprintf(fp, "\t $(RANLIB) $(UAMMlib)\n");
   fprintf(fp, "\t touch %clib.grd\n", pre);
   fprintf(fp, "clean : %cclean\n", pre);
   fprintf(fp, "%cclean:\n\t- rm -f $(objs)\n", pre);
   fprintf(fp, "killall : %ckillall\n", pre);
   fprintf(fp, "%ckillall : %cclean\n", pre, pre);
   fprintf(fp, "\t- $(ARCHIVER) d $(UAMMlib) $(objs)\n");
   fprintf(fp, "\t $(RANLIB) $(UAMMlib)\n");
   fprintf(fp, "\t- rm -f ATL_%c*.[S,c]\n", pre);

/*
 * Print out the individual rules for each needed copy function
 */
   dcomp[2] = dflags[2] = toupper(pre);
   dflags[3] = dcomp[3] = 'K';
   fprintf(fp, "\ntsth.o : tsth.c\n");
   fprintf(fp, "\t%s %s $(CDEFS) %s -c tsth.c\n\n", dcomp, dflags, styp);
   fprintf(fp, "#\n# Data copy rules\n#\n");
/*
 * Print out 2-D ablk2Cmat, cmat2ablk and ammswp targets
 */
   for (up=mnurb; up; up = up->next)
   {
      const int mu=up->mu, nu = up->nu;
      char cbe[4] = {'0', '1', 'n', 'X'};
      int ibe[4] =  {0,    1,  -1,  2};
      int i, j;

      for (i=0; i < 4; i++)
      {
         int h;
         for (h=0; h < 2; h++)
         {
            if (pre != pres[h])
               continue;
            if (i == 0)
            {
               fprintf(fp, "ATL_%cammswp%dx%d.o : ATL_%cammswp%dx%d.c\n",
                       pre, mu, nu, pre, mu, nu);
               fprintf(fp, "\t%s %s %s $(CDEFS) -c ATL_%cammswp%dx%d.c\n",
                       dcomp, dflags, styp, pre, mu, nu);
            }
            char pre=pres[h];
            char *styp=styps[h];
            for (j=0; j < 3; j++)
            {
               int k;
               char *rtn[2] = {"ablk2cmat", "cmat2ablk"};
               for (k=0; k < 2; k++)
               {
                  char rn[64];
                  sprintf(rn, "ATL_%cu%d%s_%dx%d_a%c_b%c",
                          pre, UID, rtn[k], mu, nu, al[j], cbe[i]);
                  fprintf(fp, "%s.o : %s.c\n", rn, rn);
                  fprintf(fp, "\t%s %s $(CDEFS) %s -c -DATL_%c%s=%s \\\n",
                          dcomp, dflags, styp, rn[4], rn+6+UIL, rn);
                  fprintf(fp, "          -c %s.c\n", rn);
               }
            }
         }
      }
   }
   KillUnrollList(mnurb);
/*
 * Print out 1-D copy-in routine rules
 */
   for (up=allub; up; up = up->next)
   {
      const int u = up->mu, kmaj=up->kmaj;
      int j;
      for (j=0; j < 3; j++)
      {
         int h;
         for (h=0; h < 2; h++)
         {
            if (pre != pres[h]) continue;
            char *styp=styps[h];
            char *cst[2] = {"", "C"};
            char *cdef[2] = {"", "-DConj_=1 "};
            int g;
            const int G = (pre == 'c' || pre == 'z') ? 2:1;
            for (g=0; g < G; g++)
            {
               int kk;
               int nmI, nmI0=5;
               char *sp;
               char rout[32], rt0[32];
               if (kmaj > 1)
               {
                  sprintf(rout, "ATL_%cu%drm2am_a%c_%dx%d",
                          pre, UID, al[j], kmaj, u);
                  sprintf(rt0, "ATL_%crm2am_a%c_%dx%d", pre, al[j], kmaj, u);
               }
               else
               {
                  sprintf(rout, "ATL_%cu%drm2am_a%c_%d", pre, UID, al[j], u);
                  sprintf(rt0, "ATL_%crm2am_a%c_%d", pre, al[j], u);
               }
               sp = strstr(rout, "rm2am");
               nmI = sp - rout;
               for (kk=0; kk < 4; kk++)
               {
                  if (kk == 1)
                     rt0[5] = rout[6+UIL] = 'c';
                  else if (kk == 2)
                  {
                     rt0[5] = rout[6+UIL] = 'a';
                     rt0[8] = rout[9+UIL] = 'r';
                  }
                  else if (kk == 3)
                      rt0[8] = rout[9+UIL] = 'c';

                  fprintf(fp, "%s%s.o : %s.c\n", rout, cst[g], rout);
               fprintf(fp,
                       "\t%s %s $(CDEFS) %s%s -D%s%s=%s%s -o %s%s.o -c %s.c\n",
                       dcomp, dflags, cdef[g], styp, rt0, cst[g], rout, cst[g],
                       rout, cst[g], rout);
               }
            }
         }
      }
   }
   KillUnrollList(allub);
/*
 * Print out the individual rules for each kernel compile
 */
   dflags[3] = dcomp[3] = 'M';
   fprintf(fp, "#\n#  AMM kernel rules\n#\n");
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
/*
 *    Kernels that take runtime K are only compiled once, so print rules on
 *    only the first encounter of that ID
 */
      if (FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
      {
         const int id = mmp->ID;
         for (p=mmb; p != mmp; p = p->next)
            if (p->ID == id && FLAG_IS_SET(p->flag, MMF_KRUNTIME))
               break;
         if (p != mmp)
            continue;
      }
      comp = GetKernComp(mmp, dcomp, dflags, &cflags);
/*
 *    ATL_<pre>UAMMM_<ID>_<kb>_<flg>_<mu>x<nu>x<ku>_b<X>,
 *    kerns in all 3 beta cases
 */
      for (i=0; i < 3; i++)
         PrintKernComp(fp, pre, mmp, UID, comp, cflags, styp, be[i], bes[i]);
   }
/*
 * K-cleanup needs only the kernel compile rule
 */
   fprintf(fp, "#\n#  K-cleanup rules\n#\n");
   for (mmp=ukb; mmp; mmp = mmp->next)
   {
      comp = GetKernComp(mmp, dcomp, dflags, &cflags);
      for (i=0; i < 3; i++)
         PrintKernComp(fp, pre, mmp, UID, comp, cflags, styp, be[i], bes[i]);
   }
   fclose(fp);
}

void GenKerns(char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *ukb, ATL_mmnode_t *urb,
              char *outd)
/*
 * Creates/copies all required matmul kernels into outd using specified names
 */
{
   ATL_mmnode_t *mmp, *p;
   mnur_t *mnurb=NULL, *allub, *up;
   char *ln=NULL;
   int lnlen=0, dlen;
   char al[3] = {'1', 'n', 'X'};
   int ial[3] = {1,   -1,   2};
   char pres[2];
   pres[0] = pre;
   pres[1] = (pre == 's') ? 'c' : 'z';

    GenAmmSwapFiles(pre, mmb, outd);
/*
 * Build list of all unique MU/NU combos for copy routines
 */
   mnurb = GetUniqueMNUnrolls(mmb, NULL);
   mnurb = GetUniqueMNUnrolls(urb, mnurb);
   allub = GetUniqueMUnrolls(mmb, NULL);
   allub = GetUniqueMUnrolls(urb, allub);
   allub = GetUniqueNUnrolls(mmb, allub);
   allub = GetUniqueNUnrolls(urb, allub);
   dlen = strlen(outd);
/*
 * Extract every unique block-copy routine
 */
   for (up=mnurb; up; up = up->next)
   {
      const int mu=up->mu, nu=up->nu;
      char cbe[4] = {'0', '1', 'n', 'X'};
      int ibe[4] =  {0,    1,  -1,  2};
      int i, j;
      j = 64+8 + strlen(outd);
      j = (j > 128) ? j : 128;
      if (lnlen < j)
      {
         free(ln);
         lnlen = j;
         ln = malloc(j*sizeof(char));
         assert(ln);
      }
      for (i=0; i < 4; i++)
      {
         for (j=0; j < 3; j++)
         {
            char rn[64];
            int ierr;
            int k;
            for (k=0; k < 2; k++)
            {
               int h=0;
                  if (!k)
                     sprintf(rn, "ATL_%cablk2cmat_%dx%d_a%c_b%c.c",
                             pre, mu, nu, al[j], cbe[i]);
                  else
                     sprintf(rn, "ATL_%ccmat2ablk_%dx%d_a%c_b%c.c",
                             pre, mu, nu, al[j], cbe[i]);
                  sprintf(ln,
                     "make %s pre=%c mu=%d nu=%d al=%c be=%c alpha=%d beta=%d",
                          rn, pre, mu, nu, al[j], cbe[i], ial[j], ibe[i]);
                  ierr = system(ln);
                  if (ierr)
                  {
                     fprintf(stderr, "FAILED CMND='%s'\n", ln);
                     exit(ierr);
                  }
                  sprintf(ln, "mv %s %s/ATL_%cu%d%s", rn, outd, pre, UID, rn+5);
                  ierr = system(ln);
                  if (ierr)
                  {
                     fprintf(stderr, "FAILED CMND='%s'\n", ln);
                     exit(ierr);
                  }
            }
         }
      }
   }
   KillUnrollList(mnurb);

/*
 * Routines to copy back and forth from A and B
 */
   for (up=allub; up; up = up->next)
   {
      const int u = up->mu, kmaj=up->kmaj;
      int j;
      j = 16 * 40 + (strlen(outd)<<1);
      j = (j > 90) ? j : 90;
      if (lnlen < j)
      {
         free(ln);
         lnlen = j;
         ln = malloc(j*sizeof(char));
         assert(ln);
      }
      for (j=0; j < 3; j++)
      {
         int h;
         for (h=0; h < 2; h++)
         {
            int ierr;
            if (h)
               continue;
            sprintf(ln,
                    "make ATL_%crm2am_a%c_%d.c ATL_%ccm2am_a%c_%d.c "
                    "ATL_%cam2rm_a%c_%d.c ATL_%cam2cm_a%c_%d.c "
                    "pre=%c UR=%d alpha=%d al=%c kmaj=%d",
                    pre, al[j], u, pre, al[j], u,
                    pre, al[j], u, pre, al[j], u,
                    pre, u, ial[j], al[j], kmaj);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            if (kmaj > 1)
               sprintf(ln,
                       "mv ATL_%crm2am_a%c_%d.c %s/ATL_%cu%drm2am_a%c_%dx%d.c",
                       pre, al[j], u, outd, pre, UID, al[j], kmaj, u);
            else
               sprintf(ln,
                       "mv ATL_%crm2am_a%c_%d.c %s/ATL_%cu%drm2am_a%c_%d.c",
                       pre, al[j], u, outd, pre, UID, al[j], u);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            if (kmaj > 1)
               sprintf(ln,
                       "mv ATL_%ccm2am_a%c_%d.c %s/ATL_%cu%dcm2am_a%c_%dx%d.c",
                       pre, al[j], u, outd, pre, UID, al[j], kmaj, u);
            else
               sprintf(ln,
                       "mv ATL_%ccm2am_a%c_%d.c %s/ATL_%cu%dcm2am_a%c_%d.c",
                       pre, al[j], u, outd, pre, UID, al[j], u);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            if (kmaj > 1)
               sprintf(ln,
                       "mv ATL_%cam2rm_a%c_%d.c %s/ATL_%cu%dam2rm_a%c_%dx%d.c",
                       pre, al[j], u, outd, pre, UID, al[j], kmaj, u);
            else
               sprintf(ln,
                       "mv ATL_%cam2rm_a%c_%d.c %s/ATL_%cu%dam2rm_a%c_%d.c",
                       pre, al[j], u, outd, pre, UID, al[j], u);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
            if (kmaj > 1)
               sprintf(ln,
                       "mv ATL_%cam2cm_a%c_%d.c %s/ATL_%cu%dam2cm_a%c_%dx%d.c",
                       pre, al[j], u, outd, pre, UID, al[j], kmaj, u);
            else
               sprintf(ln,
                       "mv ATL_%cam2cm_a%c_%d.c %s/ATL_%cu%dam2cm_a%c_%d.c",
                       pre, al[j], u, outd, pre, UID, al[j], u);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
         }
      }
   }
   KillUnrollList(allub);
/*
 * Copy/generate every unique file, but only once
 */
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      const int id=mmp->ID;
/*
 *    For generated files, just overwrite dups with same file, won't hurt
 */
      if (id == 0)
      {
         assert(mmp->genstr);
         assert(!system(mmp->genstr));
      }
      else  /* user-supplied files copied from AMMCASES directory */
      {
/*
 *       If this is the first time we've seen this ID, it must be copied
 */
         for (p=mmb; p != mmp && p->ID != id; p = p->next);
         if (p == mmp)
         {
            int i, ierr;
            i = strlen(mmp->rout) + dlen + 16;
            if (i > lnlen)
            {
               if (ln)
                  free(ln);
               ln = malloc(i*sizeof(char));
               assert(ln);
               lnlen = i;
            }
            sprintf(ln, "cp AMMCASES/%s %s/.", mmp->rout, outd);
            ierr = system(ln);
            if (ierr)
            {
               fprintf(stderr, "FAILED CMND='%s'\n", ln);
               exit(ierr);
            }
         }
      }
   }
/*
 * Copy/generate k-cleanup which is known to be unique
 */
   for (mmp=ukb; mmp; mmp = mmp->next)
   {
      const int id=mmp->ID;
/*
 *    For generated files, just generate it using genstr
 */
      if (id == 0)
      {
         assert(mmp->genstr);
         assert(!system(mmp->genstr));
      }
      else  /* user-supplied files copied from AMMCASES directory */
      {
         int i, ierr;
         i = strlen(mmp->rout) + dlen + 16;
         if (i > lnlen)
         {
            if (ln)
               free(ln);
            ln = malloc(i*sizeof(char));
            assert(ln);
            lnlen = i;
         }
         sprintf(ln, "cp AMMCASES/%s %s/.", mmp->rout, outd);
         ierr = system(ln);
         if (ierr)
         {
            fprintf(stderr, "FAILED CMND='%s'\n", ln);
            exit(ierr);
         }
      }
   }
   if (ln)
      free(ln);
}
void GenAllFiles(char pre, ATL_mmnode_t *mmb, ATL_mmnode_t *ukb,
                 ATL_mmnode_t *kcb, ATL_mmnode_t *rkb, ATL_mmnode_t *urb,
                 ATL_mmnode_t *sqb, ATL_mmnode_t *usb,
                 char *outd)
{
   GenHeaderFiles(pre, mmb, ukb, kcb, rkb, urb, sqb, usb, outd);
   GenMakefile(pre, mmb, ukb, rkb, urb, usb, outd);
   GenKerns(pre, mmb, ukb, urb, outd);
}


int KernelInList(ATL_mmnode_t *mmb, ATL_mmnode_t *p)
/*
 * RETURNS: 1 if p is duplicated in mmb, else 0
 */
{
   ATL_mmnode_t *mp;
   if (!p || !mmb)
      return(0);
   for (mp=mmb; mp; mp = mp->next)
      if (KernelIsSame(mp, p))
         return(1);
    return(0);
}

ATL_mmnode_t *StripNonUniqueKs(ATL_mmnode_t *ukb, ATL_mmnode_t *mmb)
/*
 * Deletes any ukb node that also appears in mmb,
 * RETURNS: possibly shortened ukb
 */
{
   ATL_mmnode_t *mp, *prev;
/*
 * Delete any repetitive nodes starting unique queue
 */
   while(ukb && KernelInList(mmb, ukb))
      ukb = KillMMNode(ukb);
   if (!ukb)
      return(ukb);
/*
 * Now, delete any internal non-unique K-cleanup
 */
   prev = ukb;
   mp = ukb->next;
   while (mp)
   {
      if (KernelInList(mmb, mp))
         mp = prev->next = KillMMNode(mp);
      else
      {
         prev = mp;
         mp = mp->next;
      }
   }
   return(ukb);
}
ATL_mmnode_t *StripExactMatchKs(ATL_mmnode_t *ukb, ATL_mmnode_t *mmb)
/*
 * Deletes any ukb node that also appears in mmb with same K-value,
 * RETURNS: possibly shortened ukb
 */
{
   ATL_mmnode_t *mp, *prev;
/*
 * Delete any repetitive nodes starting unique queue
 */
   while(ukb && ExactKernelInList(mmb, ukb))
      ukb = KillMMNode(ukb);
   if (!ukb)
      return(ukb);
/*
 * Now, delete any internal non-unique K-cleanup
 */
   prev = ukb;
   mp = ukb->next;
   while (mp)
   {
      if (ExactKernelInList(mmb, mp))
         mp = prev->next = KillMMNode(mp);
      else
      {
         prev = mp;
         mp = mp->next;
      }
   }
   return(ukb);
}
/*
 * RETURNS: mmb, with all repeated kernels removed
 */
ATL_mmnode_t *RemoveNonUniqueKernels(ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mp;
   if (!mmb || !mmb->next)
      return(mmb);
   for (mp=mmb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p=mp->next, *prev=mp;
      while (p)
      {
         if (KernelIsSame(mp, p))
            prev->next = p = KillMMNode(p);
         else
         {
            prev = p;
            p = p->next;
         }
      }
   }
   return(mmb);
}

int main(int nargs, char **args)
{
   char pre='d';
   int verb=1;
   int *nbs;
   char *outd, *ukin, *kcin, *rkin, *sqin;
   ATL_mmnode_t *mmb, *mmp, *ukb=NULL, *kcb=NULL, *rkb=NULL, *urb=NULL,
                *sqb=NULL, *usb=NULL;

   mmb = GetFlags(nargs, args, &pre, &outd, &ukin, &kcin, &rkin, &sqin);
   assert(mmb);
   if (ukin)
   {
      ukb = ReadMMFile(ukin);
      free(ukin);
/*
 *    Now, strip off any ukb node that also appears in mmb, in order to
 *    make sure we don't repeat prototypes, generation, etc
 */
      ukb = StripNonUniqueKs(ukb, mmb);
   }
   if (kcin)
   {
      kcb = ReadMMFile(kcin);
      free(kcin);
   }
   if (rkin)
   {
      rkb = ReadMMFile(rkin);
      urb = CloneMMQueue(rkb);
      urb = RemoveNonUniqueKernels(urb);
      urb = StripExactMatchKs(urb, kcb);
      urb = StripExactMatchKs(urb, mmb);
      free(rkin);
   }
   if (sqin)
   {
      sqb = ReadMMFile(sqin);
      usb = CloneMMQueue(sqb);
      usb = RemoveNonUniqueKernels(usb);
      usb = StripExactMatchKs(usb, mmb);
      usb = StripExactMatchKs(usb, kcb);
      usb = StripExactMatchKs(usb, urb);
      free(sqin);
   }
   FillInGenStrings(pre, mmb, outd);
   FillInGenStrings(pre, ukb, outd);
   FillInGenStrings(pre, kcb, outd);
   FillInGenStrings(pre, rkb, outd);
   FillInGenStrings(pre, urb, outd);
   FillInGenStrings(pre, usb, outd);
   GenAllFiles(pre, mmb, ukb, kcb, rkb, urb, sqb, usb, outd);
   KillAllMMNodes(mmb);
   KillAllMMNodes(ukb);
   KillAllMMNodes(kcb);
   KillAllMMNodes(rkb);
   exit(0);
}
