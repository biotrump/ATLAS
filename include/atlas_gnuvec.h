#ifndef ATLAS_GNUVEC
   #define ATLAS_GNUVEC 1
   #ifndef TYPE
      #if defined(SREAL) || defined(SCPLX)
         #define TYPE float
      #else
         #define TYPE double
      #endif
   #endif
   #ifdef ATL_VSX
      #define ATL_NVREG 64
      #define ATL_VLENb 16
      #if defined(SREAL) || defined (SCPLX)
         #define ATL_VLEN 4
      #else
         #define ATL_VLEN 2
      #endif
   #elif defined(ATL_AltiVec)
      #define ATL_NVREG 32
      #define ATL_VLENb 16
      #if defined(SREAL) || defined (SCPLX)
         #define ATL_VLEN 4
      #else
         #define ATL_VLEN 2
      #endif
   #elif defined(ATL_AVXZ)
      #include "immintrin.h"
      #define ATL_NVREG 32
      #define ATL_VLENb 64
      #if defined(SREAL) || defined (SCPLX)
         #define ATL_VLEN 16
         #define ATL_gvbcast(ptr_, v_) \
            v_ = _mm512_extload_ps((void*)(ptr_), _MM_UPCONV_PS_NONE, \
                                   _MM_BROADCAST_1X16, 0)
      #else
         #define ATL_VLEN 8
         #define ATL_gvbcast(ptr_, v_) \
            v_ = _mm512_extload_pd((void*)(ptr_), _MM_UPCONV_PD_NONE, \
                                   _MM_BROADCAST_1X8, 0)
      #endif
   #elif defined(ATL_AVXMAC) || defined(ATL_AVXFMA4) || defined(ATL_AVX)
      #ifdef ATL_GAS_x8664
         #define ATL_NVREG 16
      #else
         #define ATL_NVREG 8
      #endif
      #define ATL_VLENb 32
      #if defined(SREAL) || defined (SCPLX)
         #define ATL_VLEN 8
      #else
         #define ATL_VLEN 4
      #endif
      #define ATL_gvbcast(ptr_, v_) \
         (v_) = __builtin_ia32_vbroadcastsd256((void*)(ptr_));
   #elif defined(ATL_SSE3) || defined(ATL_SSE2)
      #ifdef ATL_GAS_x8664
         #define ATL_NVREG 16
      #else
         #define ATL_NVREG 8
      #endif
      #define ATL_VLENb 16
      #if defined(SREAL) || defined (SCPLX)
         #define ATL_VLEN 4
      #else
         #define ATL_VLEN 2
         #ifdef ATL_SSE3
            #define ATL_gvbcast(ptr_, v_) \
               (v_) = __builtin_ia32_loaddup(ptr_);
         #endif
      #endif
   #elif defined(SREAL) || defined(SCPLX)   /* single-only stuff */
      #ifdef ATL_AltiVec
         #define ATL_NVREG 32
         #define ATL_VLENb 16
         #define ATL_VLEN 4
      #elif defined(ATL_SSE1)
         #ifdef ATL_GAS_x8664
            #define ATL_NVREG 16
         #else
            #define ATL_NVREG 8
         #endif
         #define ATL_VLENb 16
         #define ATL_VLEN 4
      #elif defined(ATL_NONIEEE) && ATL_NONIEEE != 0
         #ifdef ATL_NEON
            #define ATL_NVREG 16
            #define ATL_VLENb 8
            #define ATL_VLEN 2
         #elif defined(ATL_3DNow)
            #define ATL_NVREG 8
            #define ATL_VLENb 16
            #define ATL_VLEN 4
         #endif
      #endif
   #endif
   #ifndef ATL_VLENb
      #define ATL_VLEN 1
      #if defined(SREAL) || defined (SCPLX)
         #define ATL_VLENb 4
      #else
         #define ATL_VLENb 8
      #endif
      #if defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)
         #define ATL_NVREG 8
      #else
         #define ATL_NVREG 32
      #endif
   #endif
   #ifndef ATL_vec_t
      #if ATL_VLEN > 1
         typedef TYPE ATL_vec_t  __attribute__ ((vector_size (ATL_VLENb)));
      #else
         #define ATL_vec_t TYPE
      #endif
   #endif
/*
 * Setup macros to multiply and divide by VLEN using shifts
 */
   #if ATL_VLEN == 1
      #define ATL_DivByVLEN(i_) (i_)
      #define ATL_MulByVLEN(i_) (i_)
   #elif ATL_VLEN == 2
      #define ATL_DivByVLEN(i_) ((i_)>>1)
      #define ATL_MulByVLEN(i_) ((i_)<<1)
   #elif ATL_VLEN == 4
      #define ATL_DivByVLEN(i_) ((i_)>>2)
      #define ATL_MulByVLEN(i_) ((i_)<<2)
   #elif ATL_VLEN == 8
      #define ATL_DivByVLEN(i_) ((i_)>>3)
      #define ATL_MulByVLEN(i_) ((i_)<<3)
   #elif ATL_VLEN == 16
      #define ATL_DivByVLEN(i_) ((i_)>>4)
      #define ATL_MulByVLEN(i_) ((i_)<<4)
   #elif ATL_VLEN == 32
      #define ATL_DivByVLEN(i_) ((i_)>>5)
      #define ATL_MulByVLEN(i_) ((i_)<<5)
   #else
      #define ATL_DivByVLEN(i_) ((i_)/ATL_VLEN)
      #define ATL_MulByVLEN(i_) ((i_)*ATL_VLEN)
   #endif
   #ifndef ATL_gvbcast
      #if ATL_VLEN == 1
         #define ATL_gvbcast(ptr_, v_) \
            { ATL_vec_t z={*(ptr_)}; v_ = z; }
      #elif ATL_VLEN == 2
         #define ATL_gvbcast(ptr_, v_) \
            { ATL_vec_t z={*(ptr_),*(ptr_)}; v_ = z; }
      #elif ATL_VLEN == 4
         #define ATL_gvbcast(ptr_, v_) \
            { ATL_vec_t z={*(ptr_),*(ptr_),*(ptr_),*(ptr_)}; v_ = z; }
      #elif ATL_VLEN == 8
         #define ATL_gvbcast(ptr_, v_) \
            v_ = {*(ptr_), *(ptr_), *(ptr_), *(ptr_), \
                  *(ptr_), *(ptr_), *(ptr_), *(ptr_)}
      #elif ATL_VLEN == 16
         #define ATL_gvbcast(ptr_, v_) \
         { \
            ATL_vec_t z_ = {*(ptr_), *(ptr_), *(ptr_), *(ptr_), \
                  *(ptr_), *(ptr_), *(ptr_), *(ptr_), \
                  *(ptr_), *(ptr_), *(ptr_), *(ptr_), \
                  *(ptr_), *(ptr_), *(ptr_), *(ptr_), \
                 } \
            v_ = z_; \
         }
      #else
         #error "Cannot create gvbcast"
      #endif
   #endif

#endif
