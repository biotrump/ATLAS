#include "atlconf.h"

#define NAFF 11
static char *AFFS[NAFF] =
   {"NONE", "SETAFFNP", "SETPROCNP", "WIN64", "WIN", "PBIND", "BINDP", "RUNON",
    "SCHED", "PLPA", "CPUSET"};
#define IA_CPUSET 10
#define IA_PLPA 9
#define IA_SCHED 8
#define IA_RUNON 7
#define IA_BINDP 6
#define IA_PBIND 5
#define IA_WIN 4
#define IA_WIN64 3
#define IA_SETPROCNP 2
#define IA_SETAFFNP 1
#define IA_NONE 0
#define IA_FIRST_SELF IA_PBIND

int RunAffProbe(int iaff, int verb, char *targ, int iproc)
{
   char ln[512], ln2[512];
   char *cmnd, *res, *frm;
   char *lib;
   int i=1;

   if (iaff == IA_WIN || iaff == IA_WIN64)
      lib = "-lkernel32";
   else
      lib = "-lpthread";
   if (targ)
   {
      i += strlen(targ);
      frm = "make IRun_%s atlrun=atlas_runX targ=%s LIBS='%s' args='%d' 2> /dev/null | fgrep SUCCESS";
   }
   else
      frm = "make IRun_%s LIBS='%s' args='%d' 2> /dev/null | fgrep SUCCESS";
   i += strlen(frm) + strlen(AFFS[iaff]) + strlen(lib) + 11;
   cmnd = malloc(sizeof(char)*i);
   assert(cmnd);
   if (targ)
      sprintf(cmnd, frm, AFFS[iaff], targ, lib, iproc);
   else
      sprintf(cmnd, frm, AFFS[iaff], lib, iproc);
   res = atlsys_1L(NULL, cmnd, verb, 0);
   if (res)
   {
      if (strstr(res, "SUCCESS"))
      {
         free(cmnd);
         free(res);
         if (verb)
            fprintf(stdout, "   %s: DETECTED!\n", AFFS[iaff]);
         return(1);
      }
   }
   if (verb > 1)
      fprintf(stdout, "   cmnd='%s' out='%s'\n", cmnd, res);
   free(cmnd);
   free(res);
   if (verb)
      fprintf(stdout, "   %s: NO.\n", AFFS[iaff]);
   return(0);
}

int GetPreferredAffinity(enum OSTYPE OS, int OMP)
/*
 * RETURNS: most likely affinity to work according to the OS reading
 */
{
   #ifdef ATL_OS_Linux
      return(OMP ? IA_SCHED : IA_SETAFFNP);
   #elif defined(ATL_OS_SunOS)
      return(IA_PBIND);
   #elif defined(ATL_OS_IRIX)
      return(IA_RUNON);
   #elif defined(ATL_OS_AIX)
      return(IA_BINDP);
   #elif defined(ATL_OS_HPUX)
      return(OMP ? IA_SCHED : IA_SETPROCNP);
   #elif defined(ATL_OS_WinNT)
      return(IA_WIN);
   #elif defined(ATL_OS_Win64)
      return(IA_WIN64);
   #else
      int iret=1;
      if (OS == OSSunOS)
         iret = IA_PBIND;
      else if (OSIsWin(OS))
         iret = IA_WIN;
      else if (OS == OSAIX)
         iret = IA_BINDP;
      else if (OS == OSHPUX)
         iret = OMP ? IA_SCHED : IA_SETPROCNP;
      else if (OSIRIX)
         iret = IA_RUNON;
      else if (OS == OSLinux && OMP)
         iret = IA_SCHED;
      return(iret);
   #endif
}

int ProbeSelfAffinity
(
   int verb,            /* verbosity */
   char *targ,          /* target machine */
   enum OSTYPE OS,      /* detected OS */
   int iaff             /* affinity used during spawn */
)
{
   int i;

   if (iaff == IA_NONE || iaff >= IA_FIRST_SELF)
      return(iaff);
   for (i=IA_FIRST_SELF; i < NAFF; i++)
      if (RunAffProbe(i, verb, targ, 0))
         return(i);
   return(IA_NONE);
}

int ProbeAffinity
(
   int verb,            /* verbosity */
   char *targ,          /* target machine */
   enum OSTYPE OS,      /* detected OS */
   int OMP              /* use OpenMP? */
)
/*
 * Searches all known affinities in AFFS for one that works to spawn to ID=0
 * RETURNS: integer indicating index in AFFS of the type of affinity
 */
{
   int ipref, iret, i;

   ipref = GetPreferredAffinity(OS, OMP);
   if (RunAffProbe(ipref, verb, targ, 0))
      return(ipref);

   for (i=1; i < NAFF; i++)
   {
/*
 *    For OpenMP, try only those methods you can use after thread startup
 */
      if (OMP && (i == IA_SETAFFNP || i == IA_SETPROCNP))
         continue;
      if (i != ipref)
      {
         if (RunAffProbe(i, verb, targ, 0))
            return(i);
      }
   }
   return(0);
}

int *ProbeAffIDs
(
   int verb,            /* verbosity */
   char *targ,          /* target machine */
   enum OSTYPE OS,      /* detected OS */
   int OMP,             /* use OpenMP? */
   int iaff             /* which entry in AFFS to search with */
)
/*
 * Finds list of IDs that can be successfully used in setting affinity
 * RETURNS: integer array, first entry is number of IDs, then IDs.
 *          If all sequential IDs work, then first entry is -maxID,
 *          and remaining entries are ignored.
 */
{
   int i, nID=0, *IDs=NULL, maxlog2, maxID, SEQ=1;

   if (!iaff)
      return(NULL);
/*
 * Find maximum power of two that works as an estimate on max assignable ID.
 * We work on the following assumptions:
 * 1. Systems that have unassignable virtual processors will have physical
 *    processors that are powers of two
 * 2. 128 is big enough to capture at least one physical processor
 * Try powers of two until at least 128.
 */
   if (verb)
      printf("FINDING LARGEST POWER OF 2 ID THAT CAN BE USED IN AFFINITY:\n");
   maxlog2 = -1;
   for (i=0; i < 8; i++)
   {
      if (RunAffProbe(iaff, 0, targ, (1<<i)))
      {
         printf("   ID=%d works.\n", (1<<i));
         maxlog2 = i;
      }
      else
         printf("   ID=%d FAILS.\n", (1<<i));
/*
 *    Windows appears to roll IDs after 32, so don't test past this point
 */
      #if defined(ATL_OS_WinNT) || defined(ATL_OS_Win9x) || \
          defined(ATL_OS_Win64)
         if (i == 5)
            break;
      #endif
   }
   if (maxlog2 == -1)
   {
      fprintf(stderr, "ONLY ID=0 WORKS FOR AFFINITY!");
      return(NULL);
   }
/*
 * If 128 worked, continue to search for even more legal assignable IDs,
 * but don't search 4096
 */
   if (maxlog2 == 7)
   {
      for (i=maxlog2+1;  i < 12; i++)
      {
         if (RunAffProbe(iaff, 0, targ, (1<<i)))
         {
            printf("   ID=%d works.\n", (1<<i));
            maxlog2 = i;
         }
         else
            break;
      }
   }
   if (verb)
      printf("LARGEST SUCCESSFUL POWER OF 2 = %d\n\n", (1<<maxlog2));
/*
 * maxlog2 is now the largest (1<<maxlog2) that we successfully can use
 * affinity on.  Note that non-power2 systems like the 3 processor phenom,
 * 6-node MIPS, etc, will still detect something (2 & 4, respectively).
 * Estimate maxID as one beneath the next maxlog2 setting in order to
 * handle these non-power2 cases.
 */
   maxID = (1 << (maxlog2+1));
   IDs = malloc(sizeof(int)*maxID);
   assert(IDs);
/*
 * Now simply try all IDs between these regions, and record the good ones
 */
   if (verb)
      printf("FINDING ALL VALID IDs BETWEEN 0 AND %d:\n", maxID-1);
   IDs[++nID] = 0;
   SEQ = 1;             /* assume all sequential IDs work */
   for (i=1; i < maxID; i++)
   {
      if (RunAffProbe(iaff, 0, targ, i))
      {
         printf("   ID=%d works.\n", i);
         IDs[++nID] = i;
      }
      else
      {
         printf("   ID=%d FAILS.\n", i);
         SEQ = 0;
      }
   }
/*
 * Indicate all IDs valid between 0 & maxID by setting nID = -maxID
 */
   IDs[0] = (SEQ) ? -IDs[nID] : nID;
   return(IDs);
}
void PrintUsage(char *name, int iarg, char *arg)
{
   fprintf(stderr, "\nERROR around arg %d (%s).\n", iarg,
           arg ? arg : "unknown");
   fprintf(stderr, "USAGE: %s [flags] where flags are:\n", name);
   fprintf(stderr, "   -v <verb> : verbosity level\n");
   fprintf(stderr, "   -O <enum OSTYPE #>  : set OS type\n");
   fprintf(stderr, "   -o <outfile>\n");
   fprintf(stderr,
           "   -t <#> : set # of threads (-1: autodect; 0: no threading)\n");
   fprintf(stderr,
           "   -tl <#> <list> : set # of threads, use list of affinity IDs\n");
   fprintf(stderr,
      "NOTE: enum #s can be found by : make xprint_enums ; ./xprint_enums\n");
   exit(iarg);
}

void GetFlags(int nargs,                /* nargs as passed into main */
              char **args,              /* args as passed into main */
              int *verb,                /* verbosity setting */
              enum OSTYPE *OS,          /* OS to assume */
              int *nthreads,           /* # of threads */
              int **tids,              /* thread affinity ID list */
              int *omp,                /* Build OpenMP version of threading? */
              int *AntThr,             /* Build Antoine's threads? */
              char **outfile,
              char **targ             /* mach to ssh to*/
             )
{
   int i, k, k0, kn, DoInt;
   char *sp, *sp0;

   *verb = 0;
   *outfile = NULL;
   *targ = NULL;

   *OS = 0;
   *verb = 0;
   *nthreads = -1;
   *tids = NULL;
   *omp = *AntThr = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 't':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *nthreads = atoi(args[i]);
         if (args[i-1][2] == 'l')
         {
            *tids = malloc(*nthreads * sizeof(int));
            assert(*tids);
            for (k=0; k < *nthreads; k++)
            {
               if (++i >= nargs)
                  PrintUsage(args[0], i, "out of arguments");
               (*tids)[k] = atoi(args[i]);
            }
         }
         break;
      case 'O':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *OS = atoi(args[i]);
         break;
      case 'v':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *verb = atoi(args[i]);
         break;
      case 'T':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *targ = args[i];
         break;
      case 'o':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *outfile = args[i];
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
}

void PrintResults
(
   char *outfile,  /* NULL or stdout, : stdout, else filename to print to */
   int iaff,       /* detected affinity entry in AFFS */
   int iaffS,      /* self-affinity */
   int *IDs        /* valid ID list */
)
{
   FILE *fpout;
   if (!outfile || !strcmp(outfile, "stdout"))
      fpout = stdout;
   else if (!strcmp(outfile, "stderr"))
      fpout = stderr;
   else
   {
      fpout = fopen(outfile, "w");
      assert(fpout);
   }
   fprintf(fpout, "/* Generated by %s */\n", __FILE__);
   fprintf(fpout, "#ifndef ATL_TAFFINITY_H\n   #define ATL_TAFFINITY_H\n\n");

   fprintf(fpout, "#define ATL_HAS_AFFINITY %d\n", iaff);
   fprintf(fpout, "#define ATL_PAFF_%s 1\n", AFFS[iaff]);
   fprintf(fpout, "#define ATL_SPAFF_%s 1\n", AFFS[iaffS]);
   if (iaff)
   {
      if (iaff == IA_SETAFFNP || iaff == IA_SETPROCNP ||
          iaff == IA_WIN64 || iaff == IA_WIN)
         fprintf(fpout,
"#define ATL_PAFF_LAUNCH 1  /* affinity can be set during thread launch */\n");
      else
         fprintf(fpout,
"#define ATL_PAFF_SELF 1  /* affinity must be set by thr after launched */\n");
      if (IDs[0] < 0)
      {
         fprintf(fpout, "#define ATL_AFF_NUMID %d\n", -IDs[0]+1);
         fprintf(fpout, "#define ATL_AFF_IDSTRIDE 1\n");
         fprintf(fpout,
                 "#define ATL_RANK_IS_PROCESSORID 1 /* good IDs [0,%d] */\n",
                 -IDs[0]);
      }
      else
      {
         int i, stride;
/*
 *       Special MIC code maximizes # of physical cores as long as possible.
 *       This minimizes cache thrashing due to competing contexts.
 *       For now, reduce the scale to the number of physical cores, simply
 *       because we can't use all that scale with the small problems that
 *       can fit on the card.
 */
         #ifdef ATL_ARCH_XeonPHI
         {
            int RP = IDs[0] - 3;
            int p = IDs[0]-4;
            int n = p>>2;

            assert(IDs[1] == 0 && IDs[2] == 1);
            fprintf(fpout, "#define ATL_AFF_NUMID %d\n", p);
            fprintf(fpout, "static int ATL_affinityIDs[%d]\n", p);
            fprintf(fpout, "   = {1");
            for (i=1; i < n; i++)
               fprintf(fpout, ", %d", (i*4+1)%RP);
            for (i=0; i < n; i++)
               fprintf(fpout, ", %d", (i*4+2)%RP);
            for (i=0; i < n; i++)
               fprintf(fpout, ", %d", (i*4+3)%RP);
            for (i=0; i < n; i++)
               fprintf(fpout, ", %d", (i*4+4)%RP);
            stride = 4;
         }
         #else
            fprintf(fpout, "#define ATL_AFF_NUMID %d\n", IDs[0]);
            fprintf(fpout, "static int ATL_affinityIDs[%d]\n", IDs[0]);
            fprintf(fpout, "   = {%d", IDs[1]);
            stride = IDs[2] - IDs[1];
            for (i=1; i < IDs[0]; i++)
            {
               fprintf(fpout, ", %d", IDs[i+1]);
               if (IDs[i+1] - IDs[i] != stride)
                  stride = 0;
            }
         #endif
         fprintf(fpout, "};\n");
         if (stride)
            fprintf(fpout, "#define ATL_IDSTRIDE %d\n", stride);
         else
            fprintf(fpout,
"#define ATL_IDSTRIDE 0 /* valid IDs not separated by constant stride */\n");
      }
   }
   else /* no affinity detected */
   {
      fprintf(fpout, "#define ATL_NOAFFINITY 1\n");
      if (fpout != stderr && fpout != stdout)
      {
         FILE *fpo;
         fpo = fopen("res/aff.h", "w");
         fprintf(fpo, "#define ATL_TAFFINITY 0\n");
         fclose(fpo);
      }
   }
   fprintf(fpout, "\n#endif /* end multiple inclusion guard */\n");
   if (fpout != stderr && fpout != stdout)
      fclose(fpout);
}

int main(int nargs, char **args)
{
   int verb, iaff, iaff2, maxT, OMP, AntThr, *IDs=NULL, i, max;
   enum OSTYPE OS;
   char *targ, *outfile;

   GetFlags(nargs, args, &verb, &OS, &maxT, &IDs, &OMP, &AntThr,
            &outfile, &targ);
   iaff = ProbeAffinity(verb, targ, OS, OMP);
   iaff2 = ProbeSelfAffinity(verb, targ, OS, iaff);
   if (IDs)
   {
      int *ids;
      if (!iaff)
      {
         fprintf(stderr,
         "You have assigned particular IDs, but affinity does not work!\n");
         exit(-1);
      }
      ids = malloc((maxT+1)*sizeof(int));
      ids[0] = maxT;
      for (i=0; i < maxT; i++)
         ids[i+1] = IDs[i];
      free(IDs);
      IDs = ids;
   }
   if (iaff && !IDs)
   {
      IDs = ProbeAffIDs(verb, targ, OS, OMP, iaff);
      assert(IDs);
   }

   printf("AFFINITY TYPES= '%s' ('%s')\n", AFFS[iaff], AFFS[iaff2]);
   if (iaff)
   {
      if (IDs[0] < 0)
         printf("All IDs between [0, %d] are valid\n", -IDs[0]);
      else
      {
         max = IDs[0] + 1;
         printf("IDs[%d]=[%d", IDs[0], IDs[1]);
         for (i=2; i < max; i++)
            printf(", %d", IDs[i]);
         printf("]\n");
      }
   }
   PrintResults(outfile, iaff, iaff2, IDs);
   free(IDs);
   return(0);
}
