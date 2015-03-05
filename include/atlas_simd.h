#ifndef ATLAS_SIMD_H
   #define  ATLAS_SIMD_H 1
/*
 * This header files contains wrappers to allow you to use SIMD vector
 * extensions in a very simplified way in a type-independent manner.
 * First argument is always destination.
 */
#if defined(ATL_AVXMAC) || defined(ATL_AVX)
   #include <immintrin.h>
   #if defined(SREAL) || defined(SCPLX)
      #if defined(ATL_VLEN) && ATL_VLEN != 8  /* VLEN OVERRIDE! */
         #if ATL_VLEN == 4
            #define ATL_VTYPE __m128
         #else
            #error "Only VLEN==4 or 8 supported for single precision AVX!"
         #endif
      #else
         #define ATL_VTYPE __m256
         #define ATL_VLEN 8
      #endif
      #define ATL_vzero(v_) v_ = _mm256_setzero_ps()
      #define ATL_vbcast(v_, p_) \
         v_ =  _mm256_broadcast_ss(p_);
      #define ATL_vuld(v_, p_) v_ = _mm256_loadu_ps(p_)
      #define ATL_vld(v_, p_) v_ = _mm256_load_ps(p_)
      #define ATL_vust(p_, v_) _mm256_storeu_ps(p_, v_)
      #define ATL_vst(p_, v_) _mm256_store_ps(p_, v_)
      #define ATL_vadd(d_, s1_, s2_) d_ =  _mm256_add_ps(s1_, s2_)
      #define ATL_vmul(d_, s1_, s2_) d_ =  _mm256_mul_ps(s1_, s2_)
      #ifdef ATL_AVXMAC
         #define ATL_vmac(d_, s1_, s2_) \
            d_ = _mm256_fmadd_ps(s1_, s2_, d_)
      #else
         #define ATL_vmac(d_, s1_, s2_) \
         { ATL_VTYPE t_; \
            t_ = _mm256_mul_ps(s1_, s2_); \
            d_ = _mm256_add_ps(t_, d_); \
         }
      #endif
   #else        /* double precision */
      #if defined(ATL_VLEN) && ATL_VLEN != 4  /* VLEN OVERRIDE! */
         #if (ATL_VLEN == 2)
            #define ATL_VTYPE __m128d
         #else
            #error "Only VLEN==2 or 4 supported for double precision AVX!"
         #endif
      #else
         #define ATL_VTYPE __m256d
         #define ATL_VLEN 4
      #endif
      #define ATL_vzero(v_) v_ = _mm256_setzero_pd()
      #define ATL_vbcast(v_, p_) v_ =  _mm256_broadcast_sd(p_)
      #define ATL_vuld(v_, p_) v_ = _mm256_loadu_pd(p_)
      #define ATL_vld(v_, p_) v_ = _mm256_load_pd(p_)
      #define ATL_vust(p_, v_) _mm256_storeu_pd(p_, v_)
      #define ATL_vst(p_, v_) _mm256_store_pd(p_, v_)
      #define ATL_vadd(d_, s1_, s2_) d_ =  _mm256_add_pd(s1_, s2_)
      #define ATL_vmul(d_, s1_, s2_) d_ =  _mm256_mul_pd(s1_, s2_)
      #ifdef ATL_AVXMAC
         #define ATL_vmac(d_, s1_, s2_) \
            d_ = _mm256_fmadd_pd(s1_, s2_, d_)
      #else
         #define ATL_vmac(d_, s1_, s2_) \
         { ATL_VTYPE t_; \
            t_ = _mm256_mul_pd(s1_, s2_); \
            d_ = _mm256_add_pd(t_, d_); \
         }
      #endif
   #endif
#elif defined(ATL_SSE2) && (defined(DREAL) || defined(DCPLX))
   #include <xmmintrin.h>
   #define ATL_VTYPE __m128d
   #if defined(ATL_VLEN) && ATL_VLEN != 2
      #error "VLEN == 2 only supported size for double precision SSE!"
   #elif !defined(ATL_VLEN)
      #define ATL_VLEN 2
   #endif
   #define ATL_vzero(v_) v_ = _mm_setzero_pd()
   #define ATL_vbcast(v_, p_) v_ =  _mm_load1_pd(p_)
   #define ATL_vuld(v_, p_) v_ = _mm_loadu_pd(p_)
   #define ATL_vld(v_, p_) v_ = _mm_load_pd(p_)
   #define ATL_vust(p_, v_) _mm_storeu_pd(p_, v_)
   #define ATL_vst(p_, v_) _mm_store_pd(p_, v_)
   #define ATL_vadd(d_, s1_, s2_) d_ =  _mm_add_pd(s1_, s2_)
   #define ATL_vmul(d_, s1_, s2_) d_ =  _mm_mul_pd(s1_, s2_)
   #define ATL_vmac(d_, s1_, s2_) \
   { ATL_VTYPE t_; \
      t_ = _mm_mul_pd(s1_, s2_); \
      d_ = _mm_add_pd(t_, d_); \
   }
#elif defined(ATL_SSE1)
   #include <xmmintrin.h>
   #define ATL_VTYPE __m128
   #if defined(ATL_VLEN) && ATL_VLEN != 4
      #error "VLEN == 4 only supported size for single precision SSE!"
   #elif !defined(ATL_VLEN)
      #define ATL_VLEN 4
   #endif
   #define ATL_vzero(v_) v_ = _mm_setzero_pd()
   #define ATL_vbcast(v_, p_) v_ =  _mm_load1_ps(p_)
   #define ATL_vuld(v_, p_) v_ = _mm_loadu_ps(p_)
   #define ATL_vld(v_, p_) v_ = _mm_load_ps(p_)
   #define ATL_vust(p_, v_) _mm_storeu_ps(p_, v_)
   #define ATL_vst(p_, v_) _mm_store_ps(p_, v_)
   #define ATL_vadd(d_, s1_, s2_) d_ =  _mm_add_ps(s1_, s2_)
   #define ATL_vmul(d_, s1_, s2_) d_ =  _mm_mul_ps(s1_, s2_)
   #define ATL_vmac(d_, s1_, s2_) \
   { ATL_VTYPE t_; \
      t_ = _mm_mul_ps(s1_, s2_); \
      d_ = _mm_add_ps(t_, d_); \
   }
#else
   #if defined(ATL_VLEN) && ATL_VLEN != 1
      #error "For systems without vector support, only ATL_VLEN=1 supported!"
   #elif !defined(ATL_VLEN)
      #define ATL_VLEN 1
   #endif
   #define ATL_VTYPE TYPE
   #define ATL_vbcast(v_, p_) v_ =  *(p_)
   #define ATL_vuld(v_, p_) v_ = *(p_)
   #define ATL_vld(v_, p_) v_ = *(p_)
   #define ATL_vust(p_, v_) *(p_) =  v_
   #define ATL_vst(p_, v_) *(p_) =  v_
   #define ATL_vadd(d_, s1_, s2_) d_ =  s1_ + s2_
   #define ATL_vmul(d_, s1_, s2_) d_ =  s1_ * s2_
   #define ATL_vmac(d_, s1_, s2_) d_ += s1_ * s2_
#endif
   #if ATL_VLEN == 1
      #define ATL_VLSH 0
      #if defined(SREAL) || defined(SCPLX)
         #define ATL_VLENb 4
      #else
         #define ATL_VLENb 8
      #endif
   #endif
   #if ATL_VLEN == 2
      #define ATL_VLSH 1
      #if defined(SREAL) || defined(SCPLX)
         #define ATL_VLENb 8
      #else
         #define ATL_VLENb 16
      #endif
   #endif
   #if ATL_VLEN == 4
      #define ATL_VLSH 2
      #if defined(SREAL) || defined(SCPLX)
         #define ATL_VLENb 16
      #else
         #define ATL_VLENb 32
      #endif
   #endif
   #if ATL_VLEN == 8
      #define ATL_VLSH 3
      #if defined(SREAL) || defined(SCPLX)
         #define ATL_VLENb 32
      #else
         #define ATL_VLENb 64
      #endif
   #endif
   #if ATL_VLEN == 16
      #define ATL_VLSH 4
      #if defined(SREAL) || defined(SCPLX)
         #define ATL_VLENb 64
      #else
         #define ATL_VLENb 128
      #endif
   #endif
   #if ATL_VLEN == 32
      #define ATL_VLSH 5
      #if defined(SREAL) || defined(SCPLX)
         #define ATL_VLENb 128
      #else
         #define ATL_VLENb 256
      #endif
   #endif

#endif  /* end multiple-inclusion guard */
