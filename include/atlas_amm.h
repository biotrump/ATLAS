#ifndef ATLAS_AMM_H
   #define ATLAS_AMM_H

#ifndef ATL_MaxMalloc_MB
   #define ATL_MaxMalloc_MB 2048UL
#endif
#ifndef ATL_MaxMalloc
   #ifdef ATL_MaxMalloc_MB
      #define ATL_MaxMalloc (ATL_MaxMalloc_MB<<20)
   #else
      #define ATL_MaxMalloc (1<<27)   /* 128 MB */
   #endif
#endif
#include "atlas_misc.h"

#ifdef TREAL
   typedef void (*cm2am_t)(const size_t, const size_t, const SCALAR,
                           const TYPE*, const size_t, TYPE*);
   typedef void (*am2cm_t)(const size_t, const size_t, const SCALAR,
                           TYPE*, const size_t, const TYPE*);
   typedef void (*ablk2cmat_t)(const size_t, const size_t, const SCALAR,
                               const TYPE*, const SCALAR, TYPE *, const size_t);
   typedef void (*cmat2ablk_t)(const size_t, const size_t, const SCALAR,
                               const TYPE*, const size_t, const SCALAR,TYPE*);
   typedef void (*ammswp_t)(ATL_CINT, TYPE*,ATL_CSZT,TYPE*);
#else
   typedef void (*cm2am_t)(const size_t, const size_t, const SCALAR,
                           const TYPE*, const size_t, TYPE*, TYPE*);
   typedef void (*am2cm_t)(const size_t, const size_t, const SCALAR,
                           TYPE*, const size_t, const TYPE*, const TYPE*);
   typedef void (*ablk2cmat_t)(const size_t, const size_t, const SCALAR,
                               const TYPE*, const TYPE*, const SCALAR,
                               TYPE *, const size_t);
   typedef void (*cmat2ablk_t)(const size_t, const size_t, const SCALAR,
                               const TYPE*, const size_t, const SCALAR,
                               TYPE*,TYPE*);
   typedef void (*ammswp_t)(ATL_CINT, TYPE*,ATL_CSZT,TYPE*,TYPE*);
#endif
typedef void (*ammkern_t)(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                          TYPE*, const TYPE*, const TYPE*, const TYPE*);
typedef struct amminfo amminfo_t;
struct amminfo
{
   cm2am_t a2blk, b2blk;
   ablk2cmat_t Cblk2cm, Cblk2cm_b1;
   cmat2ablk_t cm2Cblk;
   ammkern_t amm_b0, amm_b1, amm_bn, amm_k1_b0, amm_k1_b1, amm_k1_bn;
   unsigned short IDX, mb, nb, kb, kbmin;
   unsigned char flag, mu, nu, ku;
};

enum ATL_AMMALG {ATL_amm1b, ATL_ammrkK, ATL_NMK};

int Mjoin(PATL,GetRankKInfo)
   (amminfo_t *out, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB, ATL_CSZT M,
    ATL_CSZT N, ATL_CSZT K, const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetAmmmInfo)
   (amminfo_t *out, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB, ATL_CSZT M,
    ATL_CSZT N, ATL_CSZT K, const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,tGetAmmmInfo)
   (amminfo_t *out, const unsigned int P, enum ATLAS_TRANS TA,
    enum ATLAS_TRANS TB, ATL_CSZT M, ATL_CSZT N, ATL_CSZT K,
    const SCALAR alpha, const SCALAR beta);
ablk2cmat_t Mjoin(PATL,tGetSyammInfo)
   (amminfo_t *out, const int P, enum ATLAS_TRANS TA, ATL_CSZT N, ATL_CSZT K,
    const SCALAR alpha, const SCALAR beta);
ablk2cmat_t Mjoin(PATL,tGetSyammInfo_K)
   (amminfo_t *out, const int P, enum ATLAS_TRANS TA, ATL_CSZT N, ATL_CSZT K);


void Mjoin(PATL,ammm)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_rk2)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmREC)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmMNK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_aliased_rkK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_tN)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_IP)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_rkK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_1b)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmNMK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT,
    const SCALAR,const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,
    TYPE*,ATL_CSZT);
#endif  /* end include file guard */
