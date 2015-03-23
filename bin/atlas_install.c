/*
 *             Automatically Tuned Linear Algebra Software v3.11.32
 * Copyright (C) 1998, 2015 R. Clint Whaley
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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "atlas_fopen.h"
#include "atlas_enum.h"
#include "atlas_r1parse.h"
#include "atlas_mvparse.h"
#include "atlas_mmparse.h"
#ifndef Mmax
   #define Mmax(x, y) ( (x) > (y) ? (x) : (y) )
#endif

char *redir="2>&1 | ./xatlas_tee";
char *fmake="make -f Makefile";
int DOSEARCH=1, REGS=32;
int QUERY=0;
FILE *fpI, *fparch;

#define Mciswspace(C) ( (((C) > 8) && ((C) < 14)) || ((C) == 32) )
#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )

#include <time.h>
#include <ctype.h>
void GetDate(int *month, int *day, int *year, int *hour, int *min)
{
   time_t tv;
   struct tm *tp;

   tv = time(NULL);
   tp = localtime(&tv);
   *month = tp->tm_mon + 1;
   *day = tp->tm_mday;
   *year = tp->tm_year + 1900;
   *hour = tp->tm_hour;
   *min = tp->tm_min;
}

long GetInt(FILE *fpin, long Default, char *spc, char *expstr)
/*
 * Gets a signed integral type from fpin.  If nothing or garbage is entered,
 * Default is returned.
 */
{
   char str[64];
   long iin;
   if (expstr) fprintf(stdout, "%sEnter %s [%ld]: ", spc, expstr, Default);
   if (fgets(str, 64, fpin) == NULL) return(Default);
   if (sscanf(str, " %ld ", &iin) != 1) return(Default);
   return(iin);
}

long GetIntRange(long Default, long Min, long Max, char *spc, char *expstr)
{
   long i;
   int keepOn=0;
   do
   {
      i = GetInt(stdin, Default, spc, expstr);
      if (i > Max)
      {
         keepOn = 1;
         fprintf(stderr, "\n%ld larger than max value of %ld.  Try again.\n\n",
                 i, Max);
      }
      else if (i < Min)
      {
         keepOn = 1;
         fprintf(stderr, "\n%ld smaller than min value of %ld.  Try again.\n\n",
                 i, Min);
      }
      else keepOn = 0;
   }
   while (keepOn);
   return(i);
}

long GetIntVer(long Default, long Min, long Max, char *spc, char *expstr)
{
   long i, j;

   do
   {
      i = GetIntRange(Default, Min, Max, spc, expstr);
      fprintf(stdout, "%s   You entered: %ld\n", spc, i);
      j = GetIntRange(0, 0, 1, spc, "1 to reenter, 0 accepts");
   }
   while(j);
   return(i);
}


void GetString(FILE *fpin, char *Default, char *spc, char *expstr,
               int len, char *str0)
/*
 * Get a string of length len, not including NULL terminator; pads
 * any extra len with NULLs
 */
{
   char str[512], *sp;
   int i;

   assert(len+1 <= 512);
   if (expstr)
   {
      if (Default) fprintf(stdout, "%sEnter %s [%s]: ", spc, expstr, Default);
      else fprintf(stdout, "%sEnter %s:", spc, expstr);
   }
   sp = fgets(str, 512, fpin);
   if ( (sp == NULL) || (str[0] == '\0') || (str[0] == '\n') )
   {
      if (Default) strcpy(str0, Default);
      else str0[0] = '\0';
      return;
   }
   str[len] = '\0';
   for (i=0; str0[i] = str[i]; i++);
   if (i) i--;
   while (Mciswspace(str0[i])) i--;
   while (++i < len) str0[i] = '\0';
   str0[i] = '\0';
}

void GetStrVer(char *def, char *spc, char *expstr, int len, char *str)
{
   int i;

   do
   {
      GetString(stdin, def, spc, expstr, len, str);
      fprintf(stdout, "%sYou have entered '%s'\n", spc, str);
      i = GetIntRange(0, 0, 1, spc, "1 to reenter, 0 to accept");
   }
   while(i);
}

int IsYes(char def, char *spc, char *expstr)
{
   char ch, ln[256];
   fprintf(stdout, "%s%s [%c]: ", spc, expstr, def);
   if (fgets(ln, 256, stdin) == NULL) ch=def;
   else if (ln[0] == '\0' || ln[0] == '\n') ch=def;
   else ch = ln[0];
   return( ((ch == 'y') || (ch == 'Y')) );
}

char GetChar(char def, char *spc, char *expstr)
{
   char ch, ln[256];
   fprintf(stdout, "%s%s [%c]: ", spc, expstr, def);
   if (fgets(ln, 256, stdin) == NULL) ch=def;
   else if (ln[0] == '\0' || ln[0] == '\n') ch=def;
   else ch = ln[0];
   return(ch);
}

int FileIsThere(char *nam)
{
   FILE *fp;

   fp = fopen(nam, "r");
   if (fp == NULL) return(0);
   fclose(fp);
   return(1);
}

#include <stdarg.h>
#define ATL_UseStringVarArgs
void ATL_mprintf(int np, ...)
/*
 * Prints same message to np output file streams
 */
{
   va_list argptr;
   FILE *fp[8];
   char *form;
   int i;
   #ifdef ATL_UseStringVarArgs
      char ln[1024];
   #endif

   if (np > 0)
   {
      va_start(argptr, np);
      assert(np <= 8);
      for (i=0; i < np; i++) fp[i] = va_arg(argptr, FILE *);
      form = va_arg(argptr, char *);
      #ifdef ATL_UseStringVarArgs
         vsprintf(ln, form, argptr);
         assert(strlen(ln) < 1024);/* sanity test only, will not stop overrun */
         va_end(argptr);
         for (i=0; i < np; i++) if (fp[i]) fprintf(fp[i], "%s", ln);
      #else
         for (i=0; i < np; i++) if (fp[i]) vfprintf(fp[i], form, argptr);
         va_end(argptr);
      #endif
   }
}

int GetFirstInt(char *ln)
{
   int i, iret=0;
   for (i=0; ln[i]; i++)
   {
      if (isdigit(ln[i]))
      {
         sscanf(ln+i, "%d", &iret);
         break;
      }
   }
   return(iret);
}

long long GetFirstLong(char *ln)
{
   int i;
   long long iret=0;
   for (i=0; ln[i]; i++)
   {
      if (isdigit(ln[i]))
      {
         sscanf(ln+i, "%Ld", &iret);
         break;
      }
   }
   return(iret);
}

double GetFirstDouble(char *ln)
/*
 * Gets a double, which begins wt digit or "." (i.e., won't get e10)
 */
{
   int i;
   double dret=0;
   for (i=0; ln[i]; i++)
   {
      if (isdigit(ln[i]))
      {
         if (i > 0 && ln[i-1] == '.') i--;
         sscanf(ln+i, "%lf", &dret);
         break;
      }
   }
   return(dret);
}

int GetLastInt(char *ln)
{
   int i, iret=0;
   for (i=0; ln[i]; i++);
   if (i > 0) for (i--; i > 0 && !isdigit(ln[i]); i--);
   if (i > 0 || (i == 0 && isdigit(ln[0])))
   {
      while(isdigit(ln[i]) && i > 0) i--;
      if (!isdigit(ln[i])) i++;
      sscanf(ln+i, "%d", &iret);
   }
   return(iret);
}

long long GetLastLong(char *ln)
{
   int i;
   long iret=0;
   for (i=0; ln[i]; i++);
   if (i > 0) for (i--; i > 0 && !isdigit(ln[i]); i--);
   if (i > 0 || (i == 0 && isdigit(ln[0])))
   {
      while(isdigit(ln[i]) && i > 0) i--;
      if (!isdigit(ln[i])) i++;
      sscanf(ln+i, "%ld", &iret);
   }
   return(iret);
}

long long GetLastLongWithRound(char *ln)
/*
 * Gets a long from end of line, assuming the line may be in decimal
 * format.  If given 1.6 it will return 2.  However, if given 1e256,
 * it will return 256!
 */
{
   int i;
   long iret=0;
   for (i=0; ln[i]; i++);
   if (i > 0)
   {
      for (i--; i > 0 && !isdigit(ln[i]); i--);
      if (i > 0 && ln[i] == '.')  /* allow skip of 1 decimal point */
         for (i--; i > 0 && !isdigit(ln[i]); i--);
   }
   if (i > 0 || (i == 0 && (isdigit(ln[0]) || ln[0] == '.')))
   {
      double d;
      while(isdigit(ln[i]) && i > 0) i--;
      if (!isdigit(ln[i])) i++;
      sscanf(ln+i, "%lf", &d);
      iret = (int)(d+0.5);
   }
   return(iret);
}

void PrintBanner(char *fnam, int START, int sec, int subsec, int subsubsec);
static void ATL_Cassert0(size_t cond, char *exp, char *logfile, int line)
{
   FILE *fperr;

   if (cond) return;

   fperr = fopen("INSTALL_LOG/ERROR.LOG", "a");
   if (logfile)
      ATL_mprintf(2, stderr, fperr,
                  "ERROR %d DURING %s!!.  CHECK %s FOR DETAILS.\n",
                  line, exp, logfile);
   else ATL_mprintf(2, stderr, fperr, "ERROR %d DURING %s!!.\n", line, exp);
   if (system("make error_report") == 0)
   {
      ATL_mprintf(2, stderr, fperr,
"Error report error_<ARCH>.tgz has been created in your top-level ATLAS\n");
      ATL_mprintf(2, stderr, fperr,
         "directory.  Be sure to include this file in any help request.\n");
   }
   if (fperr) fclose(fperr);
/*
 * Ignore ifs: fucking gcc/glibc forces me to make code less maintainable
 */
   if(system("cat ../../CONFIG/error.txt >> INSTALL_LOG/ERROR.LOG"));
   if(system("cat ../../CONFIG/error.txt"));
   exit(-1);
}
#define ATL_Cassert(cond_, exp_, logfile_) \
    ATL_Cassert0((size_t)(cond_), exp_, logfile_, __LINE__)

void GetInstLogFile(char *nam, int pre, int *muladd, int *pf, int *lat,
                    int *nb, int *mu, int *nu, int *ku, int *ForceFetch,
                    int *ifetch, int *nfetch, double *mflop)
{
   char ln[128];
   FILE *fp;

   fp = fopen(nam, "r");
   if (fp == NULL) fprintf(stderr, "file %s not found!!\n\n", nam);
   assert(fp);
   assert(fgets(ln, 128, fp));
   assert(fscanf(fp, " %d %d %d %d %d %d %d %d %d %d %lf\n",
          muladd, lat, pf, nb, mu, nu, ku, ForceFetch, ifetch, nfetch, mflop)
          == 11);
   fclose(fp);
}

void PrintBanner(char *fnam, int START, int sec, int subsec, int subsubsec)
{
   int month, day, year, hour, min;
   char *sep="*******************************************************************************\n";
   char ln[256];
   FILE *fp;

   sprintf(ln, "ERROR OPENING FILE %s\n", fnam);
   fp = fopen(fnam, "a");
   if (fp == NULL) return;
   GetDate(&month, &day, &year, &hour, &min);
   fprintf(fp, "\n%s%s%s", sep, sep, sep);
   if (START)
      fprintf(fp, "*       BEGAN ATLAS3.11.3 INSTALL OF SECTION %1d-%1d-%1d ON %02d/%02d/%04d AT %02d:%02d     *\n",
              sec, subsec, subsubsec, month, day, year, hour, min);
   else
      fprintf(fp, "*      FINISHED ATLAS3.11.3 INSTALL OF SECTION %1d-%1d-%1d ON %02d/%02d/%04d AT %02d:%02d   *\n",
                      sec, subsec, subsubsec, month, day, year, hour, min);
   fprintf(fp, "%s%s%s\n\n\n", sep, sep, sep);
   fclose(fp);
}

void PrintStartStop(FILE *fp0, FILE *fp1, int nspc, int START, int sec,
                    int subsec, int subsubsec, char *stagename)
{
   int i, month, day, year, hour, min;
   char ln[512];
   FILE *fp;

   GetDate(&month, &day, &year, &hour, &min);
   for (i=0; i < nspc; i++)
       ln[i] = ' ';
   if (START)
   {
      if (stagename)
         sprintf(ln+i, "BEGIN STAGE %d-%d-%d: %s at %02d:%02d\n",
                 sec, subsec, subsubsec, stagename, hour, min);
      else
         sprintf(ln+i, "BEGIN STAGE %d-%d-%d at %02d:%02d\n",
                 sec, subsec, subsubsec, hour, min);
   }
   else
      sprintf(ln+i, "DONE  STAGE %d-%d-%d at %02d:%02d\n",
              sec, subsec, subsubsec, hour, min);
   if (START)
   {
      if (fp0)
         fprintf(fp0, "\n\n");
      if (fp1)
         fprintf(fp1, "\n\n");
   }
   if (fp0)
      fprintf(fp0, "%s", ln);
   if (fp1)
      fprintf(fp1, "%s", ln);
}

int LnIsCont(char *ln)
/*
 * RETURNS: 1 if last non-whitespace char in ln is '\', and 0 otherwise
 */
{
   int i;
   for (i=0; ln[i]; i++);
   for (i--; Mciswspace(ln[i]) && i > 0; i--);
   return(ln[i] == '\\');
}

void GetuMMRES(char pre, int ID, int *muladd, int *lat,
               int *mu, int *nu, int *ku)
/*
 * Reads the user index file to obtain the settings for user-contributed
 * kernel ID
 */
{
   char ln[1024];
   char fnam[256], auth[256];
   FILE *fp;
   int i, n, itmp, id;

   sprintf(ln, "../tune/blas/gemm/%ccases.dsc", pre);
   fp = fopen(ln, "r");
   assert(fp);
   assert(fgets(ln, 1024, fp));  /* skip comment line */
   assert(fgets(ln, 1024, fp));  /* get number of user cases line */
   assert(sscanf(ln, " %d", &n) == 1);
/*
 * Now search file for ID
 */
   for (i=0; i < n; i++)
   {
      assert(fgets(ln, 1024, fp));
      assert(sscanf(ln, " %d %d %d %d %d %d %d %d %d %d %s \"%[^\"]", &id,
                    &itmp, &itmp, &itmp, &itmp, muladd, lat, mu, nu, ku,
                    fnam, auth) == 12);
      if (LnIsCont(ln))
      {
         assert(fgets(ln, 1024, fp));
         assert(fgets(ln, 1024, fp));
      }
      if (id == ID)
         break;
   }
   assert(id == ID);
   fclose(fp);
}
void GetMMRES(char pre, int *muladd, int *lat, int *nb, int *pref,
              int *mu, int *nu, int *ku, int *ff, int *iff, int *nf,
              double *mf, int *icase, char *ufile, char *auth, double *umf)
{
   int h, i, j, k;
   char ln[256];
   FILE *fp;

   sprintf(ln, "INSTALL_LOG/%cMMRES", pre);
   fp = fopen(ln, "r");
   assert(fp);
   assert(fgets(ln, 256, fp));
   assert(fgets(ln, 256, fp));
   assert(sscanf(ln, " %d %d %d %d %d %d %d %d %d %d %lf", muladd, lat, pref,
                 nb, mu, nu, ku, ff, iff, nf, mf) == 11);
   assert(fgets(ln, 256, fp));
   if ( fgets(ln, 256, fp) )  /* user-supplied GEMM was best */
   {
      assert(fscanf(fp, " %d %d %lf \"%[^\"]\" \"%[^\"]",
                    icase, &i, umf, ufile, auth) == 5);
   }
   else
   {
      ufile[0] = auth[0] = '\0';
      *umf = 0.0;
      *icase = -1;
   }
   fclose(fp);
}

void GoToTown(int ARCHDEF, int L1DEF, int TuneLA, int NREG)
{
   const char TR[2] = {'N','T'};
   char prec[4] = {'d', 's', 'z', 'c'}, pre, upre, *typ;
   char ln[1024], tnam[256], ln2[512], ln3[512], fnam[128];
   char *mulinst, *peakstr, *peakstr2;
   int nprec=4;
   int iL1, lat, muladd, maused, latuse, lbnreg;
   int len, i, j, ip, ia, ib, ncnb, pf;
   int FoundCE=0;
   int maxreg, latus, nb, mu, nu, ku;
   int ffetch, ifetch, nfetch;
   int DefInstall=0;
   long imf;
   int idU, maU, latU, muU, nuU, kuU, il1mul, pfA;
   double mfU, mf4x1, mf4x4, mf, mfp, mmmf, mfpeak[2], l1mul;
   FILE *fp, *fpsum, *fpabr;
   ATL_mmnode_t *mmp, *mmb;

   PrintBanner("INSTALL_LOG/SUMMARY.LOG", 1, 0, 0, 0);
   fpsum = fopen("INSTALL_LOG/SUMMARY.LOG", "a");
   ATL_Cassert(fpsum, "OPENING INSTALL_LOG/SUMMARY.LOG", NULL);

   ATL_Cassert(tmpnam(tnam), "GETTING TEMPFILE", NULL);

   if (L1DEF)
   {
      ATL_Cassert(system("make IBozoL1.grd\n")==0,
                  "USING BOZO L1 DEFAULTS", NULL);
   }
   if (ARCHDEF)
      DefInstall = !system("make IArchDef.grd\n");

   PrintStartStop(stdout, fpsum, 0, 1, 1, 0, 0, "SYSTEM PROBE/AUX COMPILE");


   sprintf(ln2, "INSTALL_LOG/Stage1.log");
   PrintBanner("INSTALL_LOG/Stage1.log", 1, 1, 0, 0);
   sprintf(ln, "%s IStage1 nreg=%d %s INSTALL_LOG/Stage1.log\n",
           fmake, NREG, redir);
   ATL_Cassert(system(ln)==0, "Stage 1 install", ln2);
   fp = fopen("INSTALL_LOG/L1CacheSize", "r");

   ATL_Cassert(fp, "CACHESIZE SEARCH", ln2);
   ATL_Cassert(fscanf(fp, "%d", &i) == 1, "CACHESIZE SEARCH", ln2);
   fclose(fp);
   fprintf(stdout, "\n\n   Level 1 cache size calculated as %dKB\n", i);
   fprintf(fpsum, "   Level 1 cache size calculated as %dKB.\n\n", i);
   iL1 = i;

   for (ip=0; ip < 2; ip++)
   {
      sprintf(ln, "INSTALL_LOG/%cMULADD", prec[ip]);
      fp = fopen(ln, "r");
      ATL_Cassert(fp, "FPU PROBE", NULL);
      ATL_Cassert(fscanf(fp, " %d", &muladd)==1, "FPU PROBE", ln2);
      ATL_Cassert(fscanf(fp, " %d", &lat)==1, "FPU PROBE", ln2);
      ATL_Cassert(fscanf(fp, " %lf", &mfpeak[ip])==1, "FPU PROBE", ln2);
      ATL_Cassert(fscanf(fp, " %d", &lbnreg)==1, "FPU PROBE", ln2);
      fclose(fp);
      if (muladd) mulinst = "Combined muladd instruction";
      else mulinst = "Separate multiply and add instructions";
      ATL_mprintf(2, stdout, fpsum, "   %cFPU: %s with %d cycle pipeline.\n",
                  prec[ip], mulinst, lat);
      ATL_mprintf(2, stdout, fpsum,
                  "         Apparent number of registers : %d\n", lbnreg);
      ATL_mprintf(2, stdout, fpsum,
                  "         Register-register performance=%.2fMFLOPS\n",
                  mfpeak[ip]);
   }
   PrintStartStop(stdout, fpsum, 0, 0, 1, 0, 0, NULL);

   PrintStartStop(stdout, fpsum, 0, 1, 2, 0, 0, "TYPE-DEPENDENT TUNING");

   for (ip=0; ip < nprec; ip++)
   {
      sprintf(ln, "INSTALL_LOG/%cPerfSumm.txt", prec[ip]);
      fpabr = fopen(ln, "w");
      ATL_Cassert((fpabr != NULL), "Unable to open abbreviation file", NULL);
      pre = prec[ip];
      switch(pre)
      {
      case 's':
         len = sizeof(float);
         typ = "SREAL";
         upre = 's';
         break;
      case 'd':
         len = sizeof(double);
         typ = "DREAL";
         upre = 'd';
         break;
      case 'c':
         len = sizeof(float);
         typ = "SCPLX";
         upre = 's';
         break;
      case 'z':
         len = sizeof(double);
         typ = "DCPLX";
         upre = 'd';
         break;
      }
      #ifdef ATL_CPUMHZ
         fprintf(fpabr, "Clock_rate=%d Mhz\n", ATL_CPUMHZ);
         mfp = ATL_CPUMHZ;
         peakstr = "of detected clock rate";
         peakstr2 = "clock";
      #else
         fprintf(fpabr, "Clock_rate=0 Mhz\n");
         if (pre == 'd' || pre == 'z') mfp = mfpeak[0];
         else mfp = mfpeak[1];
         peakstr = "of apparent peak";
         peakstr2 = "peak";
      #endif
      fprintf(fpabr, "%% clock      MFLOP  ROUTINE/PROBLEM\n");
      fprintf(fpabr, "=======  =========  ===============\n");

      sprintf(ln, "TUNING PREC=\'%c\' (precision %d of %d)",
               pre, ip+1, nprec);
      PrintStartStop(stdout, fpsum, 0, 1, 2, ip+1, 0, ln);
      PrintStartStop(stdout, fpsum, 3, 1, 2, ip+1, 1, "GEMVN TUNE");
      sprintf(fnam, "INSTALL_LOG/%cMVNK.sum", pre);
      if (!FileExists(fnam))
      {
         sprintf(ln2, "INSTALL_LOG/%cMVNTUNE.LOG", pre);
         PrintBanner(ln2, 1, 2, ip+1, 1);
         sprintf(ln, "%s %s pre=%c %s %s\n", fmake, fnam, pre, redir, ln2);
         fprintf(stdout, "%s", ln);
         ATL_Cassert(system(ln)==0, "MVNTUNE", ln2);
         ATL_Cassert(FileIsThere(fnam), "MVNTUNE", ln2);
         PrintBanner(ln2, 0, 2, ip+1, 1);
      }
      {
         ATL_mvnode_t *kb, *kp;

         kb = ReadMVFile(fnam);
         ATL_Cassert(kb, "MVNTUNE", NULL);
         ATL_MVSplitContexts(kb, &kb, NULL, NULL, NULL);
         kp = ATL_LastMVNode(kb);
         ATL_mprintf(2, fpsum, stdout,
              "      gemvN : main kernel %d:%s written by %s\n", kp->ID,
                     kp->rout?kp->rout : "Generated",
                     kp->auth?kp->auth : "R. Clint Whaley");
         ATL_mprintf(2, fpsum, stdout,
                     "            mu=%d, nu=%d, using %d Cache Elements\n",
                     kp->MU, kp->NU, kp->CacheElts);
         mf = Mmax(kp->mflop[0], kp->mflop[1]);
         mf = Mmax(mf, kp->mflop[2]);
         KillAllMVNodes(kb);
         ATL_mprintf(2, fpsum, stdout,
                     "              Performance = %.2f (%5.2f of %s)\n",
                     mf, (mf / mfp)*100.0, peakstr2);
         fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "kgemvN");
      }
      PrintStartStop(stdout, fpsum, 3, 0, 2, ip+1, 1, "GEMVN TUNE");
      PrintStartStop(stdout, fpsum, 3, 1, 2, ip+1, 2, "GEMVT TUNE");
      sprintf(fnam, "INSTALL_LOG/%cMVTK.sum", pre);
      if (!FileExists(fnam))
      {
         sprintf(ln2, "INSTALL_LOG/%cMVTTUNE.LOG", pre);
         PrintBanner(ln2, 1, 2, ip+1, 2);
         sprintf(ln, "%s %s pre=%c %s %s\n", fmake, fnam, pre, redir, ln2);
         fprintf(stdout, "%s", ln);
         ATL_Cassert(system(ln)==0, "MVTTUNE", ln2);
         ATL_Cassert(FileIsThere(fnam), "MVTTUNE", ln2);
         PrintBanner(ln2, 0, 2, ip+1, 2);
      }
      {
         ATL_mvnode_t *kb, *kp;

         kb = ReadMVFile(fnam);
         ATL_Cassert(kb, "MVTTUNE", NULL);
         ATL_MVSplitContexts(kb, &kb, NULL, NULL, NULL);
         kp = ATL_LastMVNode(kb);
         ATL_mprintf(2, fpsum, stdout,
              "      gemvT : main kernel %d:%s written by %s\n", kp->ID,
                     kp->rout?kp->rout : "Generated",
                     kp->auth?kp->auth : "R. Clint Whaley");
         ATL_mprintf(2, fpsum, stdout,
                     "            mu=%d, nu=%d, using %d Cache Elements\n",
                     kp->MU, kp->NU, kp->CacheElts);
         mf = Mmax(kp->mflop[0], kp->mflop[1]);
         mf = Mmax(mf, kp->mflop[2]);
         KillAllMVNodes(kb);
         ATL_mprintf(2, fpsum, stdout,
                     "              Performance = %.2f (%5.2f of %s)\n",
                     mf, (mf / mfp)*100.0, peakstr2);
         fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "kgemvT");
      }
      PrintStartStop(stdout, fpsum, 3, 0, 2, ip+1, 2, "GEMVT TUNE");
      PrintStartStop(stdout, fpsum, 3, 1, 2, ip+1, 3, "GER  TUNE");
      sprintf(fnam, "INSTALL_LOG/%cR2K.sum", pre);
      if (!FileExists(fnam))
      {
         sprintf(ln2, "INSTALL_LOG/%cR1TUNE.LOG", pre);
         PrintBanner(ln2, 1, 2, ip+1, 7);
         sprintf(ln, "%s %s pre=%c %s %s\n", fmake, fnam, pre, redir, ln2);
         fprintf(stdout, "%s", ln);
         ATL_Cassert(system(ln)==0, "R1TUNE", ln2);
         ATL_Cassert(FileIsThere(fnam), "R1TUNE", ln2);
         PrintBanner(ln2, 0, 2, ip+1, 3);
      }
      {
         ATL_r1node_t *r1b, *r1B;
         fnam[14] = '1';
         r1b = ReadR1File(fnam);
         ATL_Cassert(r1b, "R1TUNE", NULL);
         ATL_R1SplitContexts(r1b, &r1B, NULL, NULL, NULL);
         r1b = r1B;
         r1B = ATL_LastR1Node(r1B);
         ATL_mprintf(2, fpsum, stdout,
              "      ger1 : main kernel %d:%s written by %s\n", r1B->ID,
                     r1B->rout?r1B->rout : "Generated",
                     r1B->auth?r1B->auth : "R. Clint Whaley");
         ATL_mprintf(2, fpsum, stdout,
                     "            mu=%d, nu=%d, using %d Cache Elements\n",
                     r1B->MU, r1B->NU, r1B->CacheElts);
         mf = Mmax(r1B->mflop[0], r1B->mflop[1]);
         mf = Mmax(mf, r1B->mflop[2]);
         KillAllR1Nodes(r1b);
         ATL_mprintf(2, fpsum, stdout,
                     "              Performance = %.2f (%5.2f of %s)\n",
                     mf, (mf / mfp)*100.0, peakstr2);
         fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "kger1");
         fnam[14] = '2';
         r1b = ReadR1File(fnam);
         ATL_Cassert(r1b, "R1TUNE", NULL);
         ATL_R1SplitContexts(r1b, &r1B, NULL, NULL, NULL);
         r1b = r1B;
         r1B = ATL_LastR1Node(r1B);
         ATL_mprintf(2, fpsum, stdout,
              "      ger2 : main kernel %d:%s written by %s\n", r1B->ID,
                     r1B->rout?r1B->rout : "Generated",
                     r1B->auth?r1B->auth : "R. Clint Whaley");
         ATL_mprintf(2, fpsum, stdout,
                     "            mu=%d, nu=%d, using %d Cache Elements\n",
                     r1B->MU, r1B->NU, r1B->CacheElts);
         mf = Mmax(r1B->mflop[0], r1B->mflop[1]);
         mf = Mmax(mf, r1B->mflop[2]);
         KillAllR1Nodes(r1b);
         ATL_mprintf(2, fpsum, stdout,
                     "              Performance = %.2f (%5.2f of %s)\n",
                     mf, (mf / mfp)*100.0, peakstr2);
         fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "kger2");
      }
      PrintStartStop(stdout, fpsum, 3, 0, 2, ip+1, 3, "GER  TUNE");
/*
 *    Detection of CacheEdge currently missing
 */
/*
 *    GEMM install
 */
      PrintStartStop(stdout, fpsum, 3, 1, 2, ip+1, 4, "GEMM TUNE");
      if (pre == 'd' || pre == 's')
      {
         sprintf(ln2, "INSTALL_LOG/%cAMMTUNE.LOG", pre);
         sprintf(ln, "%s AMMinstall pre=%c %s %s\n", fmake, pre, redir, ln2);
         fprintf(stdout, "%s", ln);
         ATL_Cassert(system(ln)==0, "GEMM TUNE", ln2);
         mmb = ReadMMFileWithPath(pre, "INSTALL_LOG", "eAMMRES.sum");
         ATL_Cassert(mmb, "ACCESS-MAJOR MATMUL TUNE", ln2);
         for (mmp=mmb; mmp->next; mmp = mmp->next);
         ATL_mprintf(2, fpsum, stdout,
    "\n      The best access-major matmul kernel was %s, MB=%d, NB=%d, KB=%d\n",
                     mmp->rout, mmp->mbB, mmp->nbB, mmp->kbB);
         ATL_mprintf(2, fpsum, stdout, "      written %s\n",
                     mmp->auth ? mmp->auth : "R. Clint Whaley");
         mf = mmp->mflop[0];
         ATL_mprintf(2, fpsum, stdout,
                     "      Performance: %.2fMFLOPS (%.2f of %s)\n",
                     mf, (mf / mfp)*100.0, peakstr2);
         fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "kgemm");
         sprintf(ln, "%s AMM_use pre=%c", fmake, pre);
         ATL_Cassert(system(ln)==0, "BUILDING BLOCK AMM TUNE", ln2);
         KillAllMMNodes(mmb);
      }
      else  /* complex types use real kernels, but must retune NBs */
      {
         sprintf(ln2, "INSTALL_LOG/atlas_%camm_sum.h", pre);
         if (!FileExists(ln2))
         {
            sprintf(ln2, "INSTALL_LOG/%cAMMNBTUNE.LOG", pre);
            sprintf(ln, "%s cAMMinstall pre=%c %s %s\n", fmake, pre, redir,ln2);
            fprintf(stdout, "%s", ln);
            ATL_Cassert(system(ln)==0, "COMPLEX AMM NB REFINEMENT", ln2);
         }
         sprintf(ln, "grep ATL_CAMM_APERF INSTALL_LOG/atlas_%camm_sum.h | "
                 "sed -e \"s/#define //\" -e \"s/ATL_CAMM_APERF //\""
                 ">! INSTALL_LOG/%cGEMM.MFLOP\n", pre, pre);
         mf = 1.0;
         if (system(ln))
         {
            sprintf(ln, "INSTALL_LOG/%cGEMM.MFLOP", pre);
            fp = fopen(ln, "r");
            if (fp)
            {
               if (fscanf(fp, " %le", &mf) != 1)
                  mf = 1.0;
               fclose(fp);
            }
         }
         fprintf(fpabr, "%7.1f %10.1f  %s\n", (mf/mfp)*100.0, mf, "cgemm");
      }
      fclose(fpabr);

      fprintf(fpsum, "\n");
      PrintStartStop(stdout, fpsum, 3, 0, 2, ip+1, 1, NULL);

      sprintf(ln2, "INSTALL_LOG/%cL3TUNE.LOG", pre);
      PrintBanner(ln2, 1, 2, ip+1, 5);
      PrintStartStop(stdout, fpsum, 3, 1, 2, ip+1, 4, "L3BLAS TUNE");
      if (pre == 's' || pre == 'd')
      {
         sprintf(ln, "%s INSTALL_LOG/atlas_%ctrsmXover.h pre=%c %s %s\n",
                 fmake, pre, pre, redir, ln2);
         fprintf(stdout, "%s", ln);
         ATL_Cassert(system(ln)==0, "L3BLAS TUNING", ln2);
      }
      else
      {
         sprintf(ln, "%s Il3lib pre=%c %s %s\n", fmake, pre, redir, ln2);
         fprintf(stdout, "%s", ln);
         ATL_Cassert(system(ln)==0, "L3BLAS TUNING", ln2);
      }
      sprintf(ln, "%s %ccblaslib %s %s\n", fmake, pre, redir, ln2); /* cblas */
      fprintf(stdout, "%s", ln);
      ATL_Cassert(system(ln)==0, "L3BLAS TUNING", ln2);
      PrintBanner(ln2, 0, 2, ip+1, 5);
      PrintStartStop(stdout, fpsum, 3, 0, 2, ip+1, 4, "L3BLAS TUNE");

   }
   PrintStartStop(stdout, fpsum, 0, 0, 2, 0, 0, "TYPE-DEPENDENT TUNING");
   PrintStartStop(stdout, fpsum, 0, 1, 3, 0, 0, "GENERAL LIBRARY BUILD");

   sprintf(ln2, "INSTALL_LOG/LIBBUILD.LOG");
   PrintBanner(ln2, 1, 3, 1, 1);
   sprintf(ln, "%s IBuildLibs %s %s\n", fmake, redir, ln2);
   fprintf(stdout, "%s", ln);
   ATL_Cassert(system(ln)==0, "LIBRARY BUILD", ln2);
   ATL_Cassert(FileIsThere(fnam), "LIBRARY BUILD", ln2);
   PrintBanner(ln2, 0, 3, 1, 1);
   PrintStartStop(stdout, fpsum, 0, 0, 3, 0, 0, "GENERAL LIBRARY BUILD");

   for (ip=2; ip < 4; ip++)
   {
      char pre=prec[ip];
      PrintStartStop(stdout, fpsum, 3, 1, 3, ip-1, 2,
                     "COMPLEX AMM NB REFINEMENT");
      sprintf(ln2, "INSTALL_LOG/atlas_%camm_sum.h", pre);
      if (!FileExists(ln2))
      {
         sprintf(ln2, "INSTALL_LOG/%cAMMNBTUNE.LOG", pre);
         PrintBanner(ln2, 1, 3, ip-1, 2);
         sprintf(ln, "%s cAMMinstall pre=%c %s %s\n", fmake, pre, redir, ln2);
         fprintf(stdout, "%s", ln);
         ATL_Cassert(system(ln)==0, "COMPLEX AMM NB REFINEMENT", ln2);
         PrintBanner(ln2, 0, 3, ip-1, 2);
      }
   }
   PrintStartStop(stdout, fpsum, 0, 1, 4, 0, 0, "POST-BUILD TUNING");
   sprintf(ln2, "INSTALL_LOG/POSTTUNE.LOG");
   PrintBanner(ln2, 1, 4, 1, 1);
   PrintStartStop(stdout, fpsum, 3, 1, 4, 1, 1, "TRSM TUNE");
   sprintf(ln, "%s IPostTune %s %s\n", fmake, redir, ln2);
   fprintf(stdout, "%s", ln);
   ATL_Cassert(system(ln)==0, "POST-BUILD TUNE", ln2);
   PrintStartStop(stdout, fpsum, 3, 0, 4, 1, 0, NULL);
   ATL_Cassert(FileIsThere(fnam), "POST-BUILD TUNE", ln2);
   PrintBanner(ln2, 0, 4, 1, 1);

#ifdef ATL_NCPU
   PrintStartStop(stdout, fpsum, 3, 1, 4, 2, 0, "THREADING TUNE");
   sprintf(ln2, "INSTALL_LOG/PTTUNE.LOG");
   PrintBanner(ln2, 1, 4, 2, 0);
   sprintf(ln, "%s IPTtune %s %s\n", fmake, redir, ln2);
   fprintf(stdout, "%s", ln);
   ATL_Cassert(system(ln)==0, "THREADING TUNE", ln2);
   PrintStartStop(stdout, fpsum, 3, 0, 4, 2, 0, "THREADING TUNE");

   PrintStartStop(stdout, fpsum, 3, 1, 4, 2, 1, "THREADING BUILD");
   sprintf(ln2, "INSTALL_LOG/LIBPTBUILD.LOG");
   PrintBanner(ln2, 1, 4, 2, 1);
   sprintf(ln, "%s IBuildPtlibs %s %s\n", fmake, redir, ln2);
   fprintf(stdout, "%s", ln);
   ATL_Cassert(system(ln)==0, "PTLIBRARY BUILD", ln2);
   PrintBanner(ln2, 0, 4, 2, 1);
   PrintStartStop(stdout, fpsum, 3, 0, 4, 2, 1, "THREADING BUILD");
#endif
   if (TuneLA)
   {
      PrintStartStop(stdout, fpsum, 3, 1, 4, 3, 0, "LAPACK TUNING");
      for (ip=0; ip < nprec; ip++)
      {
         pre = prec[ip];
         sprintf(ln2, "%cLAPACK TUNING", pre);
         PrintStartStop(stdout, fpsum, 6, 1, 4, 3, ip+1, ln2);
         sprintf(ln2, "INSTALL_LOG/%cLATUNE.LOG", pre);
         PrintBanner(ln2, 1, 4, 3, ip+1);
         fprintf(stdout, "%s", ln);
         sprintf(ln, "%s ILATune pre=%c %s %s\n", fmake, pre, redir, ln2);
         ATL_Cassert(system(ln)==0, "LAPACK TUNE", ln2);
         PrintStartStop(stdout, fpsum, 6, 0, 4, 3, ip+1, NULL);
         PrintBanner(ln2, 0, 4, 3, ip+1);
      }
      PrintStartStop(stdout, fpsum, 3, 0, 4, 3, 0, "LAPACK TUNING");
   }
   PrintStartStop(stdout, fpsum, 0, 0, 4, 0, 0, "POST-BUILD TUNING");
   PrintStartStop(stdout, fpsum, 0, 1, 5, 0, 0, "FINAL LIBRARY UPDATE");
   PrintStartStop(stdout, fpsum, 3, 1, 5, 1, 0, "FINAL STATIC LIBRARY UPDATE");
   sprintf(ln2, "INSTALL_LOG/LIBUPDATE.LOG");
   PrintBanner(ln2, 1, 5, 1, 1);
   #if defined(ATL_NCPU) && ATL_NCPU > 1
      sprintf(ln, "%s IBuildLibs IBuildPtlibs0 %s %s\n", fmake, redir, ln2);
   #else
      sprintf(ln, "%s IBuildLibs %s %s\n", fmake, redir, ln2);
   #endif
   fprintf(stdout, "%s", ln);
   ATL_Cassert(system(ln)==0, "STATIC LIBRARY UPDATE", ln2);
   PrintBanner(ln2, 0, 5, 1, 1);
   PrintStartStop(stdout, fpsum, 3, 0, 5, 1, 0, "FINAL STATIC LIBRARY UPDATE");
#ifdef ATL_DYLIBS
   PrintStartStop(stdout, fpsum, 3, 1, 5, 2, 0,
                  "DYNAMIC/SHARED LIBRARY UPDATE");
   sprintf(ln2, "INSTALL_LOG/LIBDYBUILD.LOG");
   PrintBanner(ln2, 1, 5, 2, 1);
   sprintf(ln, "%s IBuildDyLibs %s %s\n", fmake, redir, ln2);
   fprintf(stdout, "%s", ln);
   ATL_Cassert(system(ln)==0, "DYLIBRARY BUILD", ln2);
   PrintBanner(ln2, 0, 5, 2, 1);
   PrintStartStop(stdout, fpsum, 3, 0, 5, 2, 0, NULL);
#endif

   fprintf(stdout, "\n\n\n\n");
   fprintf(stdout, "ATLAS install complete.  Examine \n");
   fprintf(stdout, "ATLAS/bin/<arch>/INSTALL_LOG/SUMMARY.LOG for details.\n");
   fclose(fpsum);
   PrintBanner("INSTALL_LOG/SUMMARY.LOG", 0, 0, 0, 0);
}

void PrintUsage(char *nam)
{
   fprintf(stderr, "\n\nUSAGE: %s [-a <#(use archdef)> -1 <#(bozoL1)>]\n\n",
           nam);
   exit(-1);
}

void GetFlags(int nargs, char *args[], int *ARCHDEF, int *L1DEF, int *TuneLA,
              int *NREG)
{
   int i;

   *NREG = *TuneLA = *L1DEF = 0;
   *ARCHDEF = 1;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'a':
         if (++i > nargs)
            PrintUsage(args[0]);
         *ARCHDEF = atoi(args[i]);
         break;
      case 'r':
         if (++i > nargs)
            PrintUsage(args[0]);
         *NREG = atoi(args[i]);
         break;
      case 'l':
         if (++i > nargs)
            PrintUsage(args[0]);
         *TuneLA = atoi(args[i]);
         break;
      case '1':
         if (++i > nargs)
            PrintUsage(args[0]);
         *L1DEF = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
}

int main(int nargs, char *args[])
{
   int L1DEF, ARCHDEF, TuneLA, NREG;
   GetFlags(nargs, args, &ARCHDEF, &L1DEF, &TuneLA, &NREG);
   GoToTown(ARCHDEF, L1DEF, TuneLA, NREG);
   return(0);
}
