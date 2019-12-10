#ifndef HLL_H_
#define HLL_H_
//#include "common.h"
//#include "hash.h"
#include <array>
#include <climits>
#include <limits>
#include <set>
#include <cmath>
#include <atomic>
#include <random>
#include <vector>
#include <zlib.h>
#include <cstring>
#include <thread>

#include "sseutil.h"
#include "hash.h"
#include "integral.h"
#include "immintrin.h"
#include <string.h>
/** modified by liumy 
*/


//namespace sketch { namespace hll { namespace detail {
namespace Sketch { //namespace hll { namespace detail {

// Based off https://github.com/oertl/hyperloglog-sketch-estimation-paper/blob/master/c%2B%2B/cardinality_estimation.hpp
template<typename FloatType>
static constexpr FloatType gen_sigma(FloatType x) {
    if(x == 1.) return std::numeric_limits<FloatType>::infinity();
    FloatType z(x);
    for(FloatType zp(0.), y(1.); z != zp;) {
        x *= x; zp = z; z += x * y; y += y;
        if(std::isnan(z)) {
#ifndef __CUDACC__
            std::fprintf(stderr, "[W:%s:%d] Reached nan. Returning the last usable number.\n", __PRETTY_FUNCTION__, __LINE__);
#endif
            return zp;
        }
    }
    return z;
}

template<typename FloatType>
static constexpr FloatType gen_tau(FloatType x) {
    if (x == 0. || x == 1.) {
        //std::fprintf(stderr, "x is %f\n", (float)x);
        return 0.;
    }
    FloatType z(1-x), tmp(0.), y(1.), zp(x);
    while(zp != z) {
        x = std::sqrt(x);
        zp = z;
        y *= 0.5;
        tmp = (1. - x);
        z -= tmp * tmp * y;
    }
    return z / 3.;
}

//} /* detail */ } /* hll */ } /* sketch */
}  /* sketch */


//namespace sketch {
//namespace hll {
namespace Sketch {
enum EstimationMethod: uint8_t {
    ORIGINAL       = 0,
    ERTL_IMPROVED  = 1,
    ERTL_MLE       = 2
};
/**
 * 
using hash::Type;
using hash::VType;
using hash::WangHash;
using hash::MurFinHash;
*/

static constexpr const char *JESTIM_STRINGS []
{
    "ORIGINAL", "ERTL_IMPROVED", "ERTL_MLE", "ERTL_JOINT_MLE"
};
enum JointEstimationMethod: uint8_t {
    //ORIGINAL       = 0,
    //ERTL_IMPROVED  = 1, // Improved but biased method
    //ERTL_MLE       = 2, // element-wise max, followed by MLE
    ERTL_JOINT_MLE = 3  // Ertl special version
};

static const char *EST_STRS [] {
    "original",
    "ertl_improved",
    "ertl_mle",
    "ertl_joint_mle"
};

#ifdef MANUAL_CHECKS
#  ifndef VERIFY_SUM
#    define VERIFY_SUM 1
#  endif
#  ifndef VERIFY_SIMD_JOINT
#    define VERIFY_SIMD_JOINT 1
#  endif
#  ifndef VERIFY_SIMD
#    define VERIFY_SIMD 1
#  endif
#endif


using CountArrayType = std::array<uint32_t, 64>;

//namespace detail {
template<typename T>
static double ertl_ml_estimate(const T& c, unsigned p, unsigned q, double relerr=1e-2); // forward declaration
template<typename Container>
inline std::array<uint32_t, 64> sum_counts(const Container &con);
//} //namespace detail
#if VERIFY_SIMD_JOINT
/*
 *Returns the estimated number of elements:
 * [0] uniquely in h1
 * [1] uniquely in h2
 * [2] in the intersection
 * size of the union is [0] + [1] + [2]
 * size of the intersection is [2]
 */
template<typename IType, size_t N, typename=typename std::enable_if<std::is_integral<IType>::value>::type>
std::string counts2str(const std::array<IType, N> &arr) {
    std::string ret;
    for(const auto el: arr) {
        ret += std::to_string(el);
        ret += ',';
    }
    ret.pop_back();
    return ret;
}

template<typename HllType>
std::array<double, 3> ertl_joint_simple(const HllType &h1, const HllType &h2) {
    using ertl_ml_estimate;
    std::array<double, 3> ret;
    auto p = h1.p();
    auto q = h1.q();
    auto c1 = sum_counts(h1.core());
    auto c2 = sum_counts(h2.core());
    const double cAX = ertl_ml_estimate(c1, h1.p(), h1.q());
    const double cBX = ertl_ml_estimate(c2, h2.p(), h2.q());
    //std::fprintf(stderr, "cAX ml est: %lf. cBX ml els: %lf\n", cAX, cBX);
    //const double cBX = hl2.creport();
    //const double cBX = hl2.creport();
    auto tmph = h1 + h2;
    const double cABX = tmph.report();
    std::array<uint32_t, 64> countsAXBhalf{0}, countsBXAhalf{0};
    countsAXBhalf[q] = countsBXAhalf[q] = h1.m();
    std::array<uint32_t, 64> cg1{0}, cg2{0}, ceq{0};
    {
        const auto &core1(h1.core()), &core2(h2.core());
        for(uint64_t i(0); i < core1.size(); ++i) {
            switch((core1[i] > core2[i]) << 1 | (core2[i] > core1[i])) {
                case 0:
                    ++ceq[core1[i]]; break;
                case 1:
                    ++cg2[core2[i]];
                    break;
                case 2:
                    ++cg1[core1[i]];
                    break;
                default:
                    __builtin_unreachable();
            }
        }
    }
    for(unsigned _q = 0; _q < q; ++_q) {
        // Handle AXBhalf
        countsAXBhalf[_q] = cg1[_q] + ceq[_q] + cg2[_q + 1];
        assert(countsAXBhalf[q] >= countsAXBhalf[_q]);
        countsAXBhalf[q] -= countsAXBhalf[_q];

        // Handle BXAhalf
        countsBXAhalf[_q] = cg2[_q] + ceq[_q] + cg1[_q + 1];
        assert(countsBXAhalf[q] >= countsBXAhalf[_q]);
        countsBXAhalf[q] -= countsBXAhalf[_q];
    }
    double cAXBhalf = ertl_ml_estimate(countsAXBhalf, p, q - 1);
    double cBXAhalf = ertl_ml_estimate(countsBXAhalf, p, q - 1);
#if !NDEBUG
    std::fprintf(stderr, "cAXBhalf = %lf\n", cAXBhalf);
    std::fprintf(stderr, "cBXAhalf = %lf\n", cBXAhalf);
#endif
    ret[0] = cABX - cBX;
    ret[1] = cABX - cAX;
    double cX1 = (1.5 * cBX + 1.5*cAX - cBXAhalf - cAXBhalf);
    double cX2 = 2.*(cBXAhalf + cAXBhalf) - 3.*cABX;
    ret[2] = std::max(0., 0.5 * (cX1 + cX2));
#if !NDEBUG
    std::fprintf(stderr, "Halves of contribution: %lf, %lf. Initial est: %lf. Result: %lf\n", cX1, cX2, cABX, ret[2]);
#endif
    return ret;
}
#endif


//namespace detail {
    // Miscellaneous requirements.
static constexpr double LARGE_RANGE_CORRECTION_THRESHOLD = (1ull << 32) / 30.;
static constexpr double TWO_POW_32 = 1ull << 32;
static double small_range_correction_threshold(uint64_t m) {return 2.5 * m;}



template<typename CountArrType>
static double calculate_estimate(const CountArrType &counts,
                                 EstimationMethod estim, uint64_t m, uint32_t p, double alpha, double relerr=1e-2) {
    assert(estim <= 3);
#if ENABLE_COMPUTED_GOTO
    static constexpr void *arr [] {&&ORREST, &&ERTL_IMPROVED_EST, &&ERTL_MLE_EST};
    goto *arr[estim];
    ORREST: {
#else
    switch(estim) {
        case ORIGINAL: {
#endif
        assert(estim != ERTL_MLE);
        double sum = counts[0];
        for(unsigned i = 1; i < 64; ++i) if(counts[i]) sum += std::ldexp(counts[i], -i); // 64 - p because we can't have more than that many leading 0s. This is just a speed thing.
        //for(unsigned i = 1; i < 64 - p + 1; ++i) sum += std::ldexp(counts[i], -i); // 64 - p because we can't have more than that many leading 0s. This is just a speed thing.
        double value(alpha * m * m / sum);
        if(value < small_range_correction_threshold(m)) {
            if(counts[0]) {
#if !NDEBUG
                std::fprintf(stderr, "[W:%s:%d] Small value correction. Original estimate %lf. New estimate %lf.\n",
                             __PRETTY_FUNCTION__, __LINE__, value, m * std::log(static_cast<double>(m) / counts[0]));
#endif
                value = m * std::log(static_cast<double>(m) / counts[0]);
            }
        } else if(value > LARGE_RANGE_CORRECTION_THRESHOLD) {
            // Reuse sum variable to hold correction.
            // I do think I've seen worse accuracy with the large range correction, but I would need to rerun experiments to be sure.
            sum = -std::pow(2.0L, 32) * std::log1p(-std::ldexp(value, -32));
            if(!std::isnan(sum)) value = sum;
#if !NDEBUG
            else std::fprintf(stderr, "[W:%s:%d] Large range correction returned nan. Defaulting to regular calculation.\n", __PRETTY_FUNCTION__, __LINE__);
#endif
        }
        return value;
    }
#if ENABLE_COMPUTED_GOTO
    ERTL_IMPROVED_EST: {
#else
        case ERTL_IMPROVED: {
#endif
        static const double divinv = 1. / (2.L*std::log(2.L));
        double z = m * gen_tau(static_cast<double>((m-counts[64 - p + 1]))/static_cast<double>(m));
        for(unsigned i = 64-p; i; z += counts[i--], z *= 0.5); // Reuse value variable to avoid an additional allocation.
        z += m * gen_sigma(static_cast<double>(counts[0])/static_cast<double>(m));
        return m * divinv * m / z;
    }
#if ENABLE_COMPUTED_GOTO
    ERTL_MLE_EST: return ertl_ml_estimate(counts, p, 64 - p, relerr);
#else
    case ERTL_MLE: return ertl_ml_estimate(counts, p, 64 - p, relerr);
    default: __builtin_unreachable();
    }
#endif
}

template<typename CoreType>
struct parsum_data_t {
    std::atomic<uint64_t> *counts_; // Array decayed to pointer.
    const CoreType          &core_;
    const uint64_t              l_;
    const uint64_t             pb_; // Per-batch
};



union SIMDHolder {

public:

#define DEC_MAX(fn) static constexpr decltype(&fn) max_fn = &fn
#define DEC_MAX16(fn) static constexpr decltype(&fn) max_fn16 = &fn
#define DEC_MAX32(fn) static constexpr decltype(&fn) max_fn32 = &fn
#define DEC_MAX64(fn) static constexpr decltype(&fn) max_fn64 = &fn
#define DEC_GT(fn)  static constexpr decltype(&fn) gt_fn  = &fn
#define DEC_EQ(fn)  static constexpr decltype(&fn) eq_fn  = &fn

// vpcmpub has roughly the same latency as
// vpcmpOPub, but twice the throughput, so we use it instead.
// (the *epu8_mask version over the *epu8_mask)
#if HAS_AVX_512
    using SType    = __m512i;
    DEC_MAX16(_mm512_max_epu16);
    DEC_MAX32(_mm512_max_epu32);
    DEC_MAX64(_mm512_max_epu64);
#  if __AVX512BW__
    DEC_MAX(_mm512_max_epu8);
    using MaskType = __mmask64;
    static_assert(sizeof(MaskType) == sizeof(__mmask64), "Mask type should be 64-bits in size.");
    static_assert(sizeof(__mmask64) == sizeof(unsigned long long), "Mask type should be 64-bits in size.");
    DEC_GT(_mm512_cmpgt_epu8_mask);
    DEC_EQ(_mm512_cmpeq_epu8_mask);
#  else
    __m256i subs[2];
    using MaskType = SIMDHolder;
    static SIMDHolder gt_fn(__m512i a, __m512i b) {
        SIMDHolder ret;
        ret.subs[0] = _mm256_cmpgt_epi8(*(__m256i *)(&a), *(__m256i *)(&b));
        ret.subs[1] = _mm256_cmpgt_epi8(*(__m256i *)(((uint8_t *)&a) + 32), *(__m256i *)(((uint8_t *)&b) + 32));
#if !NDEBUG
        SIMDHolder ac(a), bc(b);
        for(unsigned i(0); i < sizeof(ret); ++i) {
            assert(!!ret.vals[i] == !!(ac.vals[i] > bc.vals[i]));
        }
#endif
        return ret;
    }
    static SIMDHolder max_fn(__m512i a, __m512i b) {
        SIMDHolder ret;
        ret.subs[0] = _mm256_max_epu8(*(__m256i *)(&a), *(__m256i *)(&b));
        ret.subs[1] = _mm256_max_epu8(*(__m256i *)(((uint8_t *)&a) + 32), *(__m256i *)(((uint8_t *)&b) + 32));
        return ret;
    }
    static SIMDHolder eq_fn(__m512i a, __m512i b) {
        SIMDHolder ret;
        ret.subs[0] = _mm256_cmpeq_epi8(*(__m256i *)(&a), *(__m256i *)(&b));
        ret.subs[1] = _mm256_cmpeq_epi8(*(__m256i *)(((uint8_t *)&a) + 32), *(__m256i *)(((uint8_t *)&b) + 32));
        return ret;
    }
#  endif
#elif __AVX2__
    using SType = __m256i;
    using MaskType = SIMDHolder;
    DEC_MAX(_mm256_max_epu8);
    DEC_MAX16(_mm256_max_epu16);
    DEC_MAX32(_mm256_max_epu32);
    DEC_EQ (_mm256_cmpeq_epi8);
    DEC_GT (_mm256_cmpgt_epi8);
    uint64_t sub64s[sizeof(__m256i) / sizeof(uint64_t)];
    static SIMDHolder max_fn64(__m256i a, __m256i b) {
        SIMDHolder ret;
        for(unsigned i = 0; i < sizeof(__m256i) / sizeof(uint64_t); ++i)
            ret.sub64s[i] = std::max(((uint64_t *)&a)[i], ((uint64_t *)&b)[i]);
        return ret;
    }
    //DEC_GT(_mm256_cmpgt_epu8_mask);
#elif __SSE2__
    using SType = __m128i;
    using MaskType = SIMDHolder;
    DEC_MAX(_mm_max_epu8);
    DEC_MAX16(_mm_max_epu16);
    DEC_MAX32(_mm_max_epu32);
    DEC_GT(_mm_cmpgt_epi8);
    DEC_EQ(_mm_cmpeq_epi8);
    uint64_t sub64s[sizeof(__m128i) / sizeof(uint64_t)];
    static SIMDHolder max_fn64(__m128i a, __m128i b) {
        SIMDHolder ret;
        for(unsigned i = 0; i < sizeof(__m128i) / sizeof(uint64_t); ++i)
            ret.sub64s[i] = std::max(((uint64_t *)&a)[i], ((uint64_t *)&b)[i]);
        return ret;
    }
#else
#  error("Need at least SSE2")
#endif
#undef DEC_MAX
#undef DEC_GT
#undef DEC_EQ
#undef DEC_MAX16
#undef DEC_MAX32
#undef DEC_MAX64

    SIMDHolder() {} // empty constructor
    SIMDHolder(SType val_) {val = val_;}
    operator SType &() {return val;}
    operator const SType &() const {return val;}
    static constexpr size_t nels  = sizeof(SType) / sizeof(uint8_t);
    static constexpr size_t nel16s  = sizeof(SType) / sizeof(uint16_t);
    static constexpr size_t nel32s  = sizeof(SType) / sizeof(uint32_t);
    static constexpr size_t nel64s  = sizeof(SType) / sizeof(uint64_t);
    static constexpr size_t nbits = sizeof(SType) / sizeof(uint8_t) * CHAR_BIT;
    using u8arr = uint8_t[nels];
    using u16arr = uint16_t[nels / 2];
    using u32arr = uint16_t[nels / 4];
    using u64arr = uint16_t[nels / 8];
    SType val;
    u8arr vals;
    u16arr val16s;
    u32arr val32s;
    u64arr val64s;
    template<typename T>
    void inc_counts(T &arr) const {
        //static_assert(std::is_same<std::decay_t<decltype(arr[0])>, uint64_t>::value, "Must container 64-bit integers.");
        unroller<T, 0, nels> ur;
        ur(*this, arr);
    }
    template<typename T, size_t iternum, size_t niter_left> struct unroller {
        void operator()(const SIMDHolder &ref, T &arr) const {
            ++arr[ref.vals[iternum]];
            unroller<T, iternum+1, niter_left-1>()(ref, arr);
        }
        void op16(const SIMDHolder &ref, T &arr) const {
            ++arr[ref.val16s[iternum]];
            unroller<T, iternum+1, niter_left-1>().op16(ref, arr);
        }
        void op32(const SIMDHolder &ref, T &arr) const {
            ++arr[ref.val32s[iternum]];
            unroller<T, iternum+1, niter_left-1>().op16(ref, arr);
        }
        void op64(const SIMDHolder &ref, T &arr) const {
            ++arr[ref.val64s[iternum]];
            unroller<T, iternum+1, niter_left-1>().op16(ref, arr);
        }
    };
    template<typename T, size_t iternum> struct unroller<T, iternum, 0> {
        void operator()(const SIMDHolder &ref, T &arr) const {}
        void op16(const SIMDHolder &ref, T &arr) const {}
        void op32(const SIMDHolder &ref, T &arr) const {}
        void op64(const SIMDHolder &ref, T &arr) const {}
    };
#define DEC_INC(nbits)\
    template<typename T>\
    void inc_counts##nbits (T &arr) const {\
        static_assert(std::is_integral<std::decay_t<decltype(arr[0])>>::value, "Counts must be integral.");\
        unroller<T, 0, nel##nbits##s> ur;\
        ur.op##nbits(*this, arr);\
    }
    DEC_INC(16)
    DEC_INC(32)
    DEC_INC(64)
#if 0
    template<typename T>
    void inc_counts16(T &arr) const {
        static_assert(std::is_integral<std::decay_t<decltype(arr[0])>>::value, "Counts must be integral.");
        unroller<T, 0, nel16s> ur;
        ur.op16(*this, arr);
    }
    template<typename T>
    void inc_arr32(T &arr) const {
        static_assert(std::is_integral<std::decay_t<decltype(arr[0])>>::value, "Counts must be integral.");
        unroller<T, 0, nel32s> ur;
        ur.op32(*this, arr);
    }
    template<typename T>
    void inc_arr64(T &arr) const {
        static_assert(std::is_integral<std::decay_t<decltype(arr[0])>>::value, "Counts must be integral.");
        unroller<T, 0, nel64s> ur;
        ur.op64(*this, arr);
    }
#endif
#undef DEC_INC
    static_assert(sizeof(SType) == sizeof(u8arr), "both items in the union must have the same size");
};

struct joint_unroller {
    using MType = SIMDHolder::MaskType;
    using SType = SIMDHolder::SType;
    // Woof....
#if defined(__AVX512BW__)
    static_assert(sizeof(MType) == sizeof(uint64_t), "Must be 64 bits");
#endif
    template<typename T, size_t iternum, size_t niter_left> struct ju_impl {
    inline void operator()(const SIMDHolder &ref1, const SIMDHolder &ref2, const SIMDHolder &u, T &arrh1, T &arrh2, T &arru, T &arrg1, T &arrg2, T &arreq, MType gtmask1, MType gtmask2, MType eqmask) const {
            ++arrh1[ref1.vals[iternum]];
            ++arrh2[ref2.vals[iternum]];
            ++arru[u.vals[iternum]];
#if __AVX512BW__
            arrg1[ref1.vals[iternum]] += gtmask1 & 1;
            arreq[ref1.vals[iternum]] += eqmask  & 1;
            arrg2[ref2.vals[iternum]] += gtmask2 & 1;
            gtmask1 >>= 1; gtmask2 >>= 1; eqmask  >>= 1;
            // TODO: Consider packing these into an SIMD type and shifting them as a set.
#else
            static_assert(sizeof(MType) == sizeof(SIMDHolder), "Wrong size?");
            arrg1[ref1.vals[iternum]] += gtmask1.vals[iternum] != 0;
            arreq[ref1.vals[iternum]] += eqmask.vals [iternum] != 0;
            arrg2[ref2.vals[iternum]] += gtmask2.vals[iternum] != 0;
#endif
            ju_impl<T, iternum+1, niter_left-1> ju;
            ju(ref1, ref2, u, arrh1, arrh2, arru, arrg1, arrg2, arreq, gtmask1, gtmask2, eqmask);
        }
    };
    template<typename T, size_t iternum> struct ju_impl<T, iternum, 0> {
        inline void operator()(const SIMDHolder &ref1, const SIMDHolder &ref2, const SIMDHolder &u, T &arrh1, T &arrh2, T &arru, T &arrg1, T &arrg2, T &arreq, MType gtmask1, MType gtmask2, MType eqmask) const {}
    };
    template<typename T>
    inline void operator()(const SIMDHolder &ref1, const SIMDHolder &ref2, const SIMDHolder &u, T &arrh1, T &arrh2, T &arru, T &arrg1, T &arrg2, T &arreq) const {
        ju_impl<T, 0, SIMDHolder::nels> ju;
#if __AVX512BW__
        auto g1 = SIMDHolder::gt_fn(ref1.val, ref2.val);
        auto g2 = SIMDHolder::gt_fn(ref2.val, ref1.val);
        auto eq = SIMDHolder::eq_fn(ref1.val, ref2.val);
#else
        auto g1 = SIMDHolder(SIMDHolder::gt_fn(ref1.val, ref2.val));
        auto g2 = SIMDHolder(SIMDHolder::gt_fn(ref2.val, ref1.val));
        auto eq = SIMDHolder(SIMDHolder::eq_fn(ref1.val, ref2.val));
#endif
        static_assert(std::is_same<MType, std::decay_t<decltype(g1)>>::value, "g1 should be the same time as MType");
        ju(ref1, ref2, u, arrh1, arrh2, arru, arrg1, arrg2, arreq, g1, g2, eq);
    }
    template<typename T>
    inline void sum_arrays(const SType *arr1, const SType *arr2, const SType *const arr1end, T &arrh1, T &arrh2, T &arru, T &arrg1, T &arrg2, T &arreq) const {
        SIMDHolder v1, v2, u;
        do {
            v1.val = *arr1++;
            v2.val = *arr2++;
            u.val  = SIMDHolder::max_fn(v1.val, v2.val);
            this->operator()(v1, v2, u, arrh1, arrh2, arru, arrg1, arrg2, arreq);
        } while(arr1 < arr1end);
    }
    template<typename T, typename VectorType>
    inline void sum_arrays(const VectorType &c1, const VectorType &c2, T &arrh1, T &arrh2, T &arru, T &arrg1, T &arrg2, T &arreq) const {
        assert(c1.size() == c2.size() || !std::fprintf(stderr, "Sizes: %zu, %zu\n", c1.size(), c2.size()));
        assert((c1.size() & (SIMDHolder::nels - 1)) == 0);
        sum_arrays(reinterpret_cast<const SType *>(&c1[0]), reinterpret_cast<const SType *>(&c2[0]), reinterpret_cast<const SType *>(&*c1.cend()), arrh1, arrh2, arru, arrg1, arrg2, arreq);
    }
};

template<typename T>
inline void inc_counts(T &counts, const SIMDHolder *p, const SIMDHolder *pend) {
    static_assert(std::is_integral<std::decay_t<decltype(counts[0])>>::value, "Counts must be integral.");
    SIMDHolder tmp;
    do {
        tmp = *p++;
        tmp.inc_counts(counts);
    } while(p < pend);
}

static inline std::array<uint32_t, 64> sum_counts(const SIMDHolder *p, const SIMDHolder *pend) {
    // Should add Contiguous Container requirement.
    std::array<uint32_t, 64> counts{0};
    inc_counts(counts, p, pend);
    return counts;
}

template<typename Container>
inline std::array<uint32_t, 64> sum_counts(const Container &con) {
    //static_assert(std::is_same<std::decay_t<decltype(con[0])>, uint8_t>::value, "Container must contain 8-bit unsigned integers.");
    return sum_counts(reinterpret_cast<const SIMDHolder *>(&*std::cbegin(con)), reinterpret_cast<const SIMDHolder *>(&*std::cend(con)));
}

//inline std::array<uint32_t, 64> sum_counts(const DefaultCompactVectorType &con) {
//    // TODO: add a check to make sure that it's doing it right
//    return sum_counts(reinterpret_cast<const SIMDHolder *>(con.get()), reinterpret_cast<const SIMDHolder *>(con.get() + con.bytes()));
//}

template<typename T, typename Container>
inline void inc_counts(T &counts, const Container &con) {
    //static_assert(std::is_same<std::decay_t<decltype(con[0])>, uint8_t>::value, "Container must contain 8-bit unsigned integers.");
    return inc_counts(counts, reinterpret_cast<const SIMDHolder *>(&*std::cbegin(con)), reinterpret_cast<const SIMDHolder *>(&*std::cend(con)));
}

template<typename CoreType>
void parsum_helper(void *data_, long index, int tid) {
    parsum_data_t<CoreType> &data(*reinterpret_cast<parsum_data_t<CoreType> *>(data_));
    uint64_t local_counts[64]{0};
    SIMDHolder tmp;
    const SIMDHolder *p(reinterpret_cast<const SIMDHolder *>(&data.core_[index * data.pb_])),
                     *pend(reinterpret_cast<const SIMDHolder *>(&data.core_[std::min(data.l_, (index+1) * data.pb_)]));
    do {
        tmp = *p++;
        tmp.inc_counts(local_counts);
    } while(p < pend);
    for(uint64_t i = 0; i < 64ull; ++i) data.counts_[i] += local_counts[i];
}

inline std::set<uint64_t> seeds_from_seed(uint64_t seed, size_t size) {
    std::mt19937_64 mt(seed);
    std::set<uint64_t> rset;
    while(rset.size() < size) rset.emplace(mt());
    return rset;
}
template<typename T>
static double ertl_ml_estimate(const T& c, unsigned p, unsigned q, double relerr) {
/*
    Note --
    Putting all these optimizations together finally gives the new cardinality estimation
    algorithm presented as Algorithm 8. The algorithm requires mainly only elementary
    operations. For very large cardinalities it makes sense to use the strong (46) instead
    of the weak lower bound (47) as second starting point for the secant method. The
    stronger bound is a much better approximation especially for large cardinalities, where
    the extra logarithm evaluation is amortized by savings in the number of iteration cycles.
   -Ertl paper.
TODO:  Consider adding this change to the method. This could improve our performance for other
*/
    const uint64_t m = 1ull << p;
    if (c[q+1] == m) return std::numeric_limits<double>::infinity();

    int kMin, kMax;
    for(kMin=0; c[kMin]==0; ++kMin);
    int kMinPrime = std::max(1, kMin);
    for(kMax=q+1; kMax && c[kMax]==0; --kMax);
    int kMaxPrime = std::min(static_cast<int>(q), kMax);
    double z = 0.;
    for(int k = kMaxPrime; k >= kMinPrime; z = 0.5*z + c[k--]);
    z = ldexp(z, -kMinPrime);
    unsigned cPrime = c[q+1];
    if(q) cPrime += c[kMaxPrime];
    double gprev;
    double x;
    double a = z + c[0];
    int mPrime = m - c[0];
    gprev = z + ldexp(c[q+1], -q); // Reuse gprev, setting to 0 after.
    x = gprev <= 1.5*a ? mPrime/(0.5*gprev+a): (mPrime/gprev)*std::log1p(gprev/a);
    gprev = 0;
    double deltaX = x;
    relerr /= std::sqrt(m);
    while(deltaX > x*relerr) {
        int kappaMinus1;
        frexp(x, &kappaMinus1);
        double xPrime = ldexp(x, -std::max(static_cast<int>(kMaxPrime+1), kappaMinus1+2));
        double xPrime2 = xPrime*xPrime;
        double h = xPrime - xPrime2/3 + (xPrime2*xPrime2)*(1./45. - xPrime2/472.5);
        for(int k = kappaMinus1; k >= kMaxPrime; --k) {
            double hPrime = 1. - h;
            h = (xPrime + h*hPrime)/(xPrime+hPrime);
            xPrime += xPrime;
        }
        double g = cPrime*h;
        for(int k = kMaxPrime-1; k >= kMinPrime; --k) {
            double hPrime = 1. - h;
            h = (xPrime + h*hPrime)/(xPrime+hPrime);
            xPrime += xPrime;
            g += c[k] * h;
        }
        g += x*a;
        if(gprev < g && g <= mPrime) deltaX *= (g-mPrime)/(gprev-g);
        else                         deltaX  = 0;
        x += deltaX;
        gprev = g;
    }
    return x*m;
}

template<typename HllType>
double ertl_ml_estimate(const HllType& c, double relerr=1e-2) {
    return ertl_ml_estimate(sum_counts(c.core()), c.p(), c.q(), relerr);
}

//} // namespace detail

template<typename HllType>
std::array<double, 3> ertl_joint(const HllType &h1, const HllType &h2) {
    assert(h1.m() == h2.m() || !std::fprintf(stderr, "sizes don't match! Size1: %zu. Size2: %zu\n", h1.size(), h2.size()));
    std::array<double, 3> ret;
    if(h1.get_jestim() != ERTL_JOINT_MLE) {
        ret[2] = h1.union_size(h2);
        ret[0] = h1.creport();
        ret[1] = h2.creport();
        ret[2] = ret[0] + ret[1] - ret[2];
        ret[0] -= ret[2];
        ret[1] -= ret[2];
        ret[2] = std::max(ret[2], 0.);
        return ret;
    }
//    using ertl_ml_estimate;
    auto p = h1.p();
    auto q = h1.q();
    std::array<uint32_t, 64> c1{0}, c2{0}, cu{0}, ceq{0}, cg1{0}, cg2{0};
    joint_unroller ju;
    ju.sum_arrays(h1.core(), h2.core(), c1, c2, cu, cg1, cg2, ceq);
    const double cAX = h1.get_is_ready() ? h1.creport() : ertl_ml_estimate(c1, h1.p(), h1.q());
    const double cBX = h2.get_is_ready() ? h2.creport() : ertl_ml_estimate(c2, h2.p(), h2.q());
    const double cABX = ertl_ml_estimate(cu, h1.p(), h1.q());
    // std::fprintf(stderr, "Made initials: %lf, %lf, %lf\n", cAX, cBX, cABX);
    std::array<uint32_t, 64> countsAXBhalf;
    std::array<uint32_t, 64> countsBXAhalf;
    countsAXBhalf[q] = h1.m();
    countsBXAhalf[q] = h1.m();
    for(unsigned _q = 0; _q < q; ++_q) {
        // Handle AXBhalf
        countsAXBhalf[_q] = cg1[_q] + ceq[_q] + cg2[_q + 1];
        assert(countsAXBhalf[q] >= countsAXBhalf[_q]);
        countsAXBhalf[q] -= countsAXBhalf[_q];

        // Handle BXAhalf
        countsBXAhalf[_q] = cg2[_q] + ceq[_q] + cg1[_q + 1];
        assert(countsBXAhalf[q] >= countsBXAhalf[_q]);
        countsBXAhalf[q] -= countsBXAhalf[_q];
    }
    double cAXBhalf = ertl_ml_estimate(countsAXBhalf, p, q - 1);
    double cBXAhalf = ertl_ml_estimate(countsBXAhalf, p, q - 1);
    //std::fprintf(stderr, "Made halves: %lf, %lf\n", cAXBhalf, cBXAhalf);
    ret[0] = cABX - cBX;
    ret[1] = cABX - cAX;
    double cX1 = (1.5 * cBX + 1.5*cAX - cBXAhalf - cAXBhalf);
    double cX2 = 2.*(cBXAhalf + cAXBhalf) - 3.*cABX;
    ret[2] = std::max(0., 0.5 * (cX1 + cX2));
    return ret;
}

template<typename HllType>
std::array<double, 3> ertl_joint(HllType &h1, HllType &h2) {
    if(h1.get_jestim() != ERTL_JOINT_MLE) h1.csum(), h2.csum();
    return ertl_joint(static_cast<const HllType &>(h1), static_cast<const HllType &>(h2));
}



constexpr double make_alpha(size_t m) {
    switch(m) {
        case 16: return .673;
        case 32: return .697;
        case 64: return .709;
        default: return 0.7213 / (1 + 1.079/m);
    }
}


// TODO: add a compact, 6-bit version
// For now, I think that it's preferable for thread safety,
// considering there's an intrinsic for the atomic load/store, but there would not
// be for bit-packed versions.

//using Type  = typename vec::SIMDTypes<uint64_t>::Type;
//using VType = typename vec::SIMDTypes<uint64_t>::VType;
//using Space = vec::SIMDTypes<uint64_t>;
//struct WangHash
//{
//    /* data */
//};
static constexpr auto AllocatorAlignment = sse::Alignment::
#if HAS_AVX_512
AVX512
#elif __AVX2__
AVX
#elif __SSE2__
SSE
#else
Normal
#pragma message("Note: no SIMD available, using scalar values")
#endif
     ; // TODO: extend for POWER9 ISA
template<typename ValueType>
using Allocator = sse::AlignedAllocator<ValueType, AllocatorAlignment>;



template<typename HashStruct=WangHash>
class hll {
// HyperLogLog implementation.
// To make it general, the actual point of entry is a 64-bit integer hash function.
// Therefore, you have to perform a hash function to convert various types into a suitable query.
// We could also cut our memory requirements by switching to only using 6 bits per element,
// (up to 64 leading zeros), though the gains would be relatively small
// given how memory-efficient this structure is.

// Attributes
protected:
    std::vector<uint8_t, Allocator<uint8_t>> core_;
    mutable double                          value_;
    uint32_t                                   np_;
    mutable uint8_t                 is_calculated_;
    EstimationMethod                        estim_;
    JointEstimationMethod                  jestim_;
    HashStruct                                 hf_;
public:
    using final_type = hll<HashStruct>;
    using HashType = HashStruct;
#if LZ_COUNTER
    std::array<std::atomic<uint64_t>, 64> clz_counts_; // To check for bias in insertion
#endif

    std::pair<size_t, size_t> est_memory_usage() const {
        return std::make_pair(sizeof(*this),
                              core_.size() * sizeof(core_[0]));
    }
    void reset() {
        std::fill(core_.begin(), core_.end(), uint64_t(0));
        value_ = 0;
        is_calculated_ = false;
    }
    uint64_t hash(uint64_t val) const {return hf_(val);}
    uint64_t m() const {return static_cast<uint64_t>(1) << np_;}
    double alpha()          const {return make_alpha(m());}
    double relative_error() const {return 1.03896 / std::sqrt(static_cast<double>(m()));}
    bool operator==(const hll &o) const {
        return np_ == o.np_ &&
               std::equal(core_.begin(), core_.end(), o.core_.begin());
    }
    // Constructor
    template<typename... Args>
    explicit hll(size_t np, EstimationMethod estim,
                       JointEstimationMethod jestim,
                       Args &&... args):
        core_(static_cast<uint64_t>(1) << np),
        value_(0.), np_(np), is_calculated_(0),
        estim_(estim), jestim_(jestim)
        , hf_(std::forward<Args>(args)...)
    {
#if LZ_COUNTER
        for(size_t i = 0; i < clz_counts_.size(); ++i)
            clz_counts_[i].store(uint64_t(0));
#endif
        //std::fprintf(stderr, "p = %u. q = %u. size = %zu\n", np_, q(), core_.size());
    }
    explicit hll(size_t np, HashStruct &&hs): hll(np, ERTL_MLE, (JointEstimationMethod)ERTL_JOINT_MLE, std::move(hs)) {}
    explicit hll(size_t np, EstimationMethod estim=ERTL_MLE): hll(np, estim, (JointEstimationMethod)ERTL_JOINT_MLE) {}
    explicit hll(): hll(0, EstimationMethod::ERTL_MLE, JointEstimationMethod::ERTL_JOINT_MLE) {}
    template<typename... Args>
    hll(const char *path, Args &&... args): hf_(std::forward<Args>(args)...) {read(path);}
    template<typename... Args>
    hll(const std::string &path, Args &&... args): hll(path.data(), std::forward<Args>(args)...) {}
    template<typename... Args>
    hll(gzFile fp, Args &&... args): hll(0, ERTL_MLE, ERTL_JOINT_MLE, std::forward<Args>(args)...) {this->read(fp);}

    // Call sum to recalculate if you have changed contents.
    void sum() const {
        const auto counts(sum_counts(core_)); // std::array<uint32_t, 64>
        value_ = calculate_estimate(counts, estim_, m(), np_, alpha());
        is_calculated_ = 1;
    }
    void csum() const {if(!is_calculated_) sum();}

    // Returns cardinality estimate. Sums if not calculated yet.
    double creport() const {
        csum();
        return value_;
    }
    auto finalize() const {return *this;}
    auto finalize() {
        auto ret(std::move(*this));
        this->free();
        return ret;
    }
    double report() noexcept {
        csum();
        return creport();
    }
    //TODO:Added by liumy to show sketch for testing. 
    void showSketch(){
        int vecSize = core_.size();
        for(int i=0; i<vecSize; ++i)
            fprintf( stderr,"%u, ",  core_[i] );
        fprintf(stderr,"\n");
    }

    double distance(const hll &h2) const {return jaccard_index(h2);}
    double cardinality_estimate() const { return creport();}
    double cardinality_estimate() noexcept { return report();}

    // Returns error estimate
    double cest_err() const {
        if(!is_calculated_) throw std::runtime_error("Result must be calculated in order to report.");
        return relative_error() * creport();
    }
    double est_err()  noexcept {
        return cest_err();
    }
    // Returns string representation
    std::string to_string() const {
        std::string params(std::string("p:") + std::to_string(np_) + '|' + EST_STRS[estim_] + ";");
        return (params + (is_calculated_ ? std::to_string(creport()) + ", +- " + std::to_string(cest_err())
                                         : desc_string()));
    }
    // Descriptive string.
    std::string desc_string() const {
        char buf[256];
        std::sprintf(buf, "Size: %u. nb: %llu. error: %lf. Is calculated: %s. value: %lf. Estimation method: %s\n",
                     np_, static_cast<long long unsigned int>(m()), relative_error(), is_calculated_ ? "true": "false", value_, EST_STRS[estim_]);
        return buf;
    }

    inline void add(uint64_t hashval) {
#ifndef NOT_THREADSAFE
        for(const uint32_t index(hashval >> q()), lzt(clz(((hashval << 1)|1) << (np_ - 1)) + 1);
            core_[index] < lzt;
            __sync_bool_compare_and_swap(core_.data() + index, core_[index], lzt));
#else
        const uint32_t index(hashval >> q());
        const uint8_t lzt(clz(((hashval << 1)|1) << (np_ - 1)) + 1);
        core_[index] = std::max(core_[index], lzt);
#endif
#if LZ_COUNTER
        ++clz_counts_[clz(((hashval << 1)|1) << (np_ - 1)) + 1];
#endif
    }

    inline void addh(uint64_t element) {
        element = hf_(element);
        add(element);
    }
    inline void addh(const std::string &element) {
#ifdef ENABLE_CLHASH
        CONST_IF(std::is_same<HashStruct, clhasher>::value) {
            add(hf_(element));
        } else {
#endif
            add(std::hash<std::string>{}(element));
#ifdef ENABLE_CLHASH
        }
#endif
    }
    void update(const char* seq) const {
        //reverse&complenment
        const uint64_t LENGTH = strlen(seq);
        char* seqRev;
        seqRev = new char[LENGTH];
		char table[4] = {'T','G','A','C'};
		for ( uint64_t i = 0; i < LENGTH; i++ )
		{
			char base = seq[i];
			base >>= 1;
			base &= 0x03;
			seqRev[LENGTH - i - 1] = table[base];
		}
        //sequence -> kmer
        //fprintf(stderr, "seqRev = %s \n", seqRev);
        const int KMERLEN = 32;
        for(uint64_t i=0; i<LENGTH-KMERLEN; ++i) {
			char kmer[KMERLEN+1];
			char kmer_fwd[KMERLEN+1];
			char kmer_rev[KMERLEN+1];
			memcpy(kmer_fwd, seq+i, KMERLEN);
			memcpy(kmer_rev, seqRev+LENGTH-i-KMERLEN, KMERLEN);
			kmer_fwd[KMERLEN] = '\0';
			kmer_rev[KMERLEN] = '\0';
            if(memcmp(kmer_fwd, kmer_rev, KMERLEN) <= 0) {
        		//fprintf(stderr, "kmer_fwd = %s \n", kmer_fwd);
				this->addh(kmer_fwd);
			} else {
        		//fprintf(stderr, "kmer_rev = %s \n", kmer_rev);
				this->addh(kmer_rev);
			}

        }
    }

//    inline void addh(VType element) {
//        element = hf_(element.simd_);
//        add(element);
//    }
//    inline void add(VType element) {
//        element.for_each([&](uint64_t &val) {add(val);});
//    }
    template<typename T, typename Hasher=std::hash<T>>
    inline void adds(const T element, const Hasher &hasher) {
        static_assert(std::is_same<std::decay_t<decltype(hasher(element))>, uint64_t>::value, "Must return 64-bit hash");
        add(hasher(element));
    }
#ifdef ENABLE_CLHASH
    template<typename Hasher=clhasher>
    inline void adds(const char *s, size_t len, const Hasher &hasher) {
        static_assert(std::is_same<std::decay_t<decltype(hasher(s, len))>, uint64_t>::value, "Must return 64-bit hash");
        add(hasher(s, len));
    }
#endif
    void parsum(int nthreads=-1, size_t pb=4096) {
        if(nthreads < 0) nthreads = nthreads > 0 ? nthreads: std::thread::hardware_concurrency();
        std::atomic<uint64_t> acounts[64];
        std::memset(acounts, 0, sizeof acounts);
        parsum_data_t<decltype(core_)> data{acounts, core_, m(), pb};
        const uint64_t nr(core_.size() / pb + (core_.size() % pb != 0));
        kt_for(nthreads, parsum_helper<decltype(core_)>, &data, nr);
        uint64_t counts[64];
        std::memcpy(counts, acounts, sizeof(counts));
        value_ = calculate_estimate(counts, estim_, m(), np_, alpha());
        is_calculated_ = 1;
    }
    ssize_t printf(std::FILE *fp) const {
        ssize_t ret = std::fputc('[', fp) > 0;
        for(size_t i = 0; i < core_.size() - 1; ++i)
            ret += std::fprintf(fp, "%d, ", int(core_[i]));
        return ret += std::fprintf(fp, "%d]", int(core_.back()));
    }
    std::string sprintf() const {
        std::fprintf(stderr, "Core size: %zu. to string: %s\n", core_.size(), to_string().data());
        std::string ret = "[";
        for(size_t i = 0; i < core_.size() - 1; ++i)
            ret += std::to_string(core_[i]), ret += ", ";
        ret += std::to_string(core_.back());
        ret += ']';
        return ret;
    }
    hll<HashStruct> compress(size_t new_np) const {
        // See Algorithm 3 in https://arxiv.org/abs/1702.01284
        // This is not very optimized.
        // I might later add support for doubling, c/o https://research.neustar.biz/2013/04/30/doubling-the-size-of-an-hll-dynamically-extra-bits/
        if(new_np == np_) return hll(*this);
        if(new_np > np_) throw std::runtime_error(std::string("Can't compress to a larger size. Current: ") + std::to_string(np_) + ". Requested new size: " + std::to_string(new_np));
        hll<HashStruct> ret(new_np, get_estim(), get_jestim());
        size_t ratio = static_cast<size_t>(1) << (np_ - new_np);
        size_t new_size = 1ull << new_np;
        size_t b = 0;
        for(size_t i(0); i < new_size; ++i) {
            size_t j(0);
            while(j < ratio && core_[j + b] == 0) ++j;
            if(j != ratio)
                ret.core_[i] = std::min(ret.q() + 1, j ? clz(j)+1: core_[b]);
            // Otherwise left at 0
            b += ratio;
        }
        return ret;
    }
    // Reset.
    void clear() {
        if(core_.size() > (1u << 16)) {
            std::memset(core_.data(), 0, core_.size() * sizeof(core_[0]));
        }// else if(__builtin_expect(core_.size() > Space::COUNT, 1)) {
          //  for(VType v1 = Space::set1(0), *p1(reinterpret_cast<VType *>(&core_[0])), *p2(reinterpret_cast<VType *>(&core_[core_.size()])); p1 < p2; *p1++ = v1);
        //} 
        else std::fill(core_.begin(), core_.end(), static_cast<uint8_t>(0));
        value_ = is_calculated_ = 0;
    }
    hll(hll&&o): value_(0), np_(0), is_calculated_(0), estim_(ERTL_MLE), jestim_(static_cast<JointEstimationMethod>(ERTL_MLE)) {
        std::swap_ranges(reinterpret_cast<uint8_t *>(this),
                         reinterpret_cast<uint8_t *>(this) + sizeof(*this),
                         reinterpret_cast<uint8_t *>(std::addressof(o)));
    }
    hll(const hll &other): core_(other.core_), value_(other.value_), np_(other.np_), is_calculated_(other.is_calculated_),
        estim_(other.estim_), jestim_(other.jestim_), hf_(other.hf_)
    {
#if LZ_COUNTER
        for(size_t i = 0; i < clz_counts_.size(); ++i)
            clz_counts_[i].store(other.clz_counts_[i].load());
#endif
    }
    hll& operator=(const hll &other) {
        // Explicitly define to make sure we don't do unnecessary reallocation.
        if(core_.size() != other.core_.size()) core_.resize(other.core_.size());
        std::memcpy(core_.data(), other.core_.data(), core_.size()); // TODO: consider SIMD copy
        np_ = other.np_;
        value_ = other.value_;
        is_calculated_ = other.is_calculated_;
        estim_ = other.estim_;
        return *this;
    }
    hll& operator=(hll&&) = default;
    hll clone() const {
        return hll(np_, estim_, jestim_);
    }

    hll &operator+=(const hll &other) {
        if(other.np_ != np_) {
            char buf[256];
            std::sprintf(buf, "For operator +=: np_ (%u) != other.np_ (%u)\n", np_, other.np_);
            throw std::runtime_error(buf);
        }
        unsigned i;
#if HAS_AVX_512 || __AVX2__ || __SSE2__
        if(m() >= sizeof(Type)) {
#if HAS_AVX_512 && __AVX512BW__
            __m512i *els(reinterpret_cast<__m512i *>(core_.data()));
            const __m512i *oels(reinterpret_cast<const __m512i *>(other.core_.data()));
            for(i = 0; i < m() >> 6; ++i) els[i] = _mm512_max_epu8(els[i], oels[i]); // mm512_max_epu8 is available on with AVX512BW :(
#elif __AVX2__
            __m256i *els(reinterpret_cast<__m256i *>(core_.data()));
            const __m256i *oels(reinterpret_cast<const __m256i *>(other.core_.data()));
            for(i = 0; i < m() * sizeof(uint8_t) / sizeof(__m256i); ++i) {
                assert(reinterpret_cast<const char *>(&els[i]) < reinterpret_cast<const char *>(&core_[core_.size()]));
                els[i] = _mm256_max_epu8(els[i], oels[i]);
            }
#else // __SSE2__
            __m128i *els(reinterpret_cast<__m128i *>(core_.data()));
            const __m128i *oels(reinterpret_cast<const __m128i *>(other.core_.data()));
            for(i = 0; i < m() >> 4; ++i) els[i] = _mm_max_epu8(els[i], oels[i]);
#endif /* #if (HAS_AVX_512 && __AVX512BW__) || __AVX2__ || true */

            if(m() < sizeof(Type)) for(;i < m(); ++i) core_[i] = std::max(core_[i], other.core_[i]);
        } else {
#endif /* #if HAS_AVX_512 || __AVX2__ || __SSE2__ */
            uint64_t *els(reinterpret_cast<uint64_t *>(core_.data()));
            const uint64_t *oels(reinterpret_cast<const uint64_t *>(other.core_.data()));
            while(els < reinterpret_cast<const uint64_t *>(core_.data() + core_.size()))
                *els = std::max(*els, *oels), ++els, ++oels;
#if HAS_AVX_512 || __AVX2__ || __SSE2__
        }
#endif
        not_ready();
        return *this;
    }

    // Clears, allows reuse with different np.
    void resize(size_t new_size) {
        if(new_size & (new_size - 1)) new_size = roundup(new_size);
        clear();
        core_.resize(new_size);
        np_ = ilog2(new_size);
    }
    EstimationMethod get_estim()       const {return  estim_;}
    JointEstimationMethod get_jestim() const {return jestim_;}
    void set_estim(EstimationMethod val) {
        estim_ = std::max(val, ERTL_MLE);
    }
    void set_jestim(JointEstimationMethod val) {
        jestim_ = val;
    }
    void set_jestim(uint16_t val) {set_jestim(static_cast<JointEstimationMethod>(val));}
    void set_estim(uint16_t val)  {estim_  = static_cast<EstimationMethod>(val);}
    // Getter for is_calculated_
    bool get_is_ready() const {return is_calculated_;}
    void not_ready() {is_calculated_ = false;}
    void set_is_ready() {is_calculated_ = true;}
    bool may_contain(uint64_t hashval) const {
        // This returns false positives, but never a false negative.
        return core_[hashval >> q()] >= clz(hashval << np_) + 1;
    }

    bool within_bounds(uint64_t actual_size) const {
        return std::abs(actual_size - creport()) < relative_error() * actual_size;
    }

    bool within_bounds(uint64_t actual_size) {
        return std::abs(actual_size - report()) < est_err();
    }
    const auto &core()    const {return core_;}
    const uint8_t *data() const {return core_.data();}

    uint32_t p() const {return np_;}
    uint32_t q() const {return (sizeof(uint64_t) * CHAR_BIT) - np_;}
    void free() {
        decltype(core_) tmp{};
        std::swap(core_, tmp);
    }
    void write(FILE *fp) const {write(fileno(fp));}
    void write(gzFile fp) const {
#define CW(fp, src, len) do {if(gzwrite(fp, src, len) == 0) throw std::runtime_error("Error writing to file.");} while(0)
        uint32_t bf[]{is_calculated_, estim_, jestim_, 1};
        CW(fp, bf, sizeof(bf));
        CW(fp, &np_, sizeof(np_));
        CW(fp, &value_, sizeof(value_));
        CW(fp, core_.data(), core_.size() * sizeof(core_[0]));
#undef CW
    }
    void write(const char *path, bool write_gz=true) const {
        if(write_gz) {
            gzFile fp(gzopen(path, "wb"));
            if(!fp) throw ZlibError(Z_ERRNO, std::string("Could not open file at ") + path);
            write(fp);
            gzclose(fp);
        } else {
            std::FILE *fp(std::fopen(path, "wb"));
            if(fp == nullptr) throw std::runtime_error(std::string("Could not open file at ") + path);
            write(fileno(fp));
            std::fclose(fp);
        }
    }
    void write(const std::string &path, bool write_gz=false) const {write(path.data(), write_gz);}
    void read(gzFile fp) {
#define CR(fp, dst, len) \
    do {\
        if(static_cast<uint64_t>(gzread(fp, dst, len)) != len) {\
            char buf[512];\
            throw std::runtime_error(std::string(buf, buf + std::sprintf(buf, "[E:%s:%d:%s] Error reading from file\n", __FILE__, __LINE__, __PRETTY_FUNCTION__))); \
        }\
    } while(0)
        uint32_t bf[4];
        CR(fp, bf, sizeof(bf));
        is_calculated_ = bf[0];
        estim_  = static_cast<EstimationMethod>(bf[1]);
        jestim_ = static_cast<JointEstimationMethod>(bf[2]);
        CR(fp, &np_, sizeof(np_));
        CR(fp, &value_, sizeof(value_));
        core_.resize(m());
        CR(fp, core_.data(), core_.size());
#undef CR
    }
    void read(const char *path) {
        gzFile fp(gzopen(path, "rb"));
        if(fp == nullptr) throw std::runtime_error(std::string("Could not open file at ") + path);
        read(fp);
        gzclose(fp);
    }
    void read(const std::string &path) {
        read(path.data());
    }
    void write(int fileno) const {
        uint32_t bf[]{is_calculated_, estim_, jestim_, 137};
#define CHWR(fn, obj, sz) if(__builtin_expect(::write(fn, (obj), (sz)) != ssize_t(sz), 0)) throw std::runtime_error(std::string("Failed to write to disk in ") + __PRETTY_FUNCTION__)
        CHWR(fileno, bf, sizeof(bf));
        CHWR(fileno, &np_, sizeof(np_));
        CHWR(fileno, &value_, sizeof(value_));
        CHWR(fileno, core_.data(), core_.size());
#undef CHWR
    }
    void read(int fileno) {
        uint32_t bf[4];
#define CHRE(fn, obj, sz) if(__builtin_expect(::read(fn, (obj), (sz)) != ssize_t(sz), 0)) throw std::runtime_error(std::string("Failed to read from fd in ") + __PRETTY_FUNCTION__)
        CHRE(fileno, bf, sizeof(bf));
        is_calculated_ = bf[0];
        estim_         = static_cast<EstimationMethod>(bf[1]);
        jestim_        = static_cast<JointEstimationMethod>(bf[2]);
#if VERBOSE_AF
        if(bf[3] != 137) {
            std::fprintf(stderr, "Warning: old sketches. Still binary compatible, but FYI.\n");
        }
#endif
        CHRE(fileno, &np_, sizeof(np_));
        CHRE(fileno, &value_, sizeof(value_));
        core_.resize(m());
        CHRE(fileno, core_.data(), core_.size());
#undef CHRE
    }
    hll operator+(const hll &other) const {
        if(other.p() != p())
            throw std::runtime_error(std::string("p (") + std::to_string(p()) + " != other.p (" + std::to_string(other.p()));
        hll ret(*this);
        ret += other;
        return ret;
    }
    double union_size(const hll &other) const {
        if(jestim_ != JointEstimationMethod::ERTL_JOINT_MLE) {
            assert(m() == other.m());
            std::array<uint32_t, 64> counts{0};
            // We can do this because we use an aligned allocator.
            // We also have found that wider vectors than SSE2 don't matter
            const __m128i *p1(reinterpret_cast<const __m128i *>(data())), *p2(reinterpret_cast<const __m128i *>(other.data()));
            const __m128i *const pe(reinterpret_cast<const __m128i *>(&core_[core_.size()]));
            for(__m128i tmp;p1 < pe;) {
                tmp = _mm_max_epu8(*p1++, *p2++);
                for(size_t i = 0; i < sizeof(tmp);++counts[reinterpret_cast<uint8_t *>(&tmp)[i++]]);
            }
            return calculate_estimate(counts, get_estim(), m(), p(), alpha());
        }
        //std::fprintf(stderr, "jestim is ERTL_JOINT_MLE: %s\n", JESTIM_STRINGS[jestim_]);
        const auto full_counts = ertl_joint(*this, other);
        return full_counts[0] + full_counts[1] + full_counts[2];
    }
    // Jaccard index, but returning a bool to indicate whether it was less than expected error for the cardinality/sketch size
    std::pair<double, bool> bjaccard_index(hll &h2) {
        if(jestim_ != JointEstimationMethod::ERTL_JOINT_MLE) csum(), h2.csum();
        return const_cast<hll &>(*this).bjaccard_index(const_cast<const hll &>(h2));
    }
    std::pair<double, bool> bjaccard_index(const hll &h2) const {
        if(jestim_ == JointEstimationMethod::ERTL_JOINT_MLE) {
            auto full_cmps = ertl_joint(*this, h2);
            auto ret = full_cmps[2] / (full_cmps[0] + full_cmps[1] + full_cmps[2]);
            return std::make_pair(ret, ret > relative_error());
        }
        const double us = union_size(h2);
        const double ret = std::max(0., creport() + h2.creport() - us) / us;
        return std::make_pair(ret, ret > relative_error());
    }
    double jaccard_index(hll &h2) {
        if(jestim_ != JointEstimationMethod::ERTL_JOINT_MLE) csum(), h2.csum();
        return const_cast<hll &>(*this).jaccard_index(const_cast<const hll &>(h2));
    }
    double containment_index(const hll &h2) const {
        auto fsr = full_set_comparison(h2);
        return fsr[2] / (fsr[2] + fsr[0]);
    }
    std::array<double, 3> full_set_comparison(const hll &h2) const {
        if(jestim_ == JointEstimationMethod::ERTL_JOINT_MLE) {
            return ertl_joint(*this, h2);
        }
        const double us = union_size(h2), mys = creport(), os = h2.creport(),
                     is = std::max(mys + os - us, 0.),
                     my_only = std::max(mys - is, 0.), o_only = std::max(os - is, 0.);
        return std::array<double, 3>{my_only, o_only, is};
    }
    double jaccard_index(const hll &h2) const {
        if(jestim_ == JointEstimationMethod::ERTL_JOINT_MLE) {
            auto full_cmps = ertl_joint(*this, h2);
            const auto ret = full_cmps[2] / (full_cmps[0] + full_cmps[1] + full_cmps[2]);
            return ret;
        }
        const double us = union_size(h2);
        const double ret = (creport() + h2.creport() - us) / us;
        return std::max(0., ret);
    }
    size_t size() const {return size_t(m());}
    static constexpr unsigned min_size() {
        return ilog2(sizeof(SIMDHolder));
    }
#if LZ_COUNTER
    ~hll() {
        std::string tmp;
        for(const auto &val: clz_counts_) tmp += std::to_string(val), tmp += ',';
        tmp.pop_back();
        std::fprintf(stderr, "counts: %s\n", tmp.data());
    }
#endif
};


// Returns the size of the set intersection
template<typename HllType>
inline double intersection_size(HllType &first, HllType &other) noexcept {
    first.csum(), other.csum();
    return intersection_size(static_cast<const HllType &>(first), static_cast<const HllType &>(other));
}

template<typename HllType> inline std::pair<double, bool> bjaccard_index(const HllType &h1, const HllType &h2) {return h1.bjaccard_index(h2);}
template<typename HllType> inline std::pair<double, bool> bjaccard_index(HllType &h1, HllType &h2) {return h1.bjaccard_index(h2);}

// Returns a HyperLogLog union
template<typename HllType>
static inline double union_size(const HllType &h1, const HllType &h2) {return h1.union_size(h2);}

template<typename HllType>
static inline double intersection_size(const HllType &h1, const HllType &h2) {
    return std::max(0., h1.creport() + h2.creport() - union_size(h1, h2));
}


//} // namespace hll


//} // namespace sketch
} // namespace Sketch

#endif // #ifndef HLL_H_
