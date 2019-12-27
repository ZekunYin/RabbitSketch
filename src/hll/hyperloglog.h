/***********************************
 *Author: liumy@mail.sdu.edu.cn
 *Date: 12/12/2019
 *Reference: Dashing
 *TODO:
 *optimization:
 *	kmer
 *	hash
 *	clz
 *  K->C
 *  union
 *function:
 *	estimation
 *	distance
 ***********************************/

#ifndef _HYPERLOGLOG_H_
#define _HYPERLOGLOG_H_

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
#include <string.h>
#include <cassert>



namespace Sketch {

	/*************clz******************************/
	// Thomas Wang hash
	// Original site down, available at https://naml.us/blog/tag/thomas-wang
	// This is our core 64-bit hash.
	// It has a 1-1 mapping from any one 64-bit integer to another
	// and can be inverted with irving_inv_hash.
	static inline std::uint64_t wang_hash(std::uint64_t key) noexcept {
		key = (~key) + (key << 21); // key = (key << 21) - key - 1;
		key = key ^ (key >> 24);
		key = (key + (key << 3)) + (key << 8); // key * 265
		key = key ^ (key >> 14);
		key = (key + (key << 2)) + (key << 4); // key * 21
		key = key ^ (key >> 28);
		key = key + (key << 31);
		return key;
	}

	static inline std::uint64_t roundup64(std::size_t x) noexcept {
		--x;
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16;
		x |= x >> 32;
		return ++x;
	}

#define clztbl(x, arg) do {\
	switch(arg) {\
		case 0:                         x += 4; break;\
		case 1:                         x += 3; break;\
		case 2: case 3:                 x += 2; break;\
		case 4: case 5: case 6: case 7: x += 1; break;\
	}} while(0)

	constexpr inline int clz_manual( std::uint32_t x )
	{
		int n(0);
		if ((x & 0xFFFF0000) == 0) {n  = 16; x <<= 16;}
		if ((x & 0xFF000000) == 0) {n +=  8; x <<=  8;}
		if ((x & 0xF0000000) == 0) {n +=  4; x <<=  4;}
		clztbl(n, x >> (32 - 4));
		return n;
	}

	// Overload
	constexpr inline int clz_manual( std::uint64_t x )
	{
		int n(0);
		if ((x & 0xFFFFFFFF00000000ull) == 0) {n  = 32; x <<= 32;}
		if ((x & 0xFFFF000000000000ull) == 0) {n += 16; x <<= 16;}
		if ((x & 0xFF00000000000000ull) == 0) {n +=  8; x <<=  8;}
		if ((x & 0xF000000000000000ull) == 0) {n +=  4; x <<=  4;}
		clztbl(n, x >> (64 - 4));
		return n;
	}

	// clz wrappers. Apparently, __builtin_clzll is undefined for values of 0.
	// For our hash function, there is only 1 64-bit integer value which causes this problem.
	// I'd expect that this is acceptable. And on Haswell+, this value is the correct value.
#if __GNUC__ || __clang__
#ifdef AVOID_CLZ_UNDEF
	constexpr inline unsigned clz(unsigned long long x) {
		return x ? __builtin_clzll(x) : sizeof(x) * CHAR_BIT;
	}

	constexpr inline unsigned clz(unsigned long x) {
		return x ? __builtin_clzl(x) : sizeof(x) * CHAR_BIT;
	}
#else
	constexpr inline unsigned clz(unsigned long long x) {
		return __builtin_clzll(x);
	}

	constexpr inline unsigned clz(unsigned long x) {
		return __builtin_clzl(x);
	}
#endif
#else
#pragma message("Using manual clz")
#define clz(x) clz_manual(x)
	// https://en.wikipedia.org/wiki/Find_first_set#CLZ
	// Modified for constexpr, added 64-bit overload.
#endif

	static_assert(clz(0x0000FFFFFFFFFFFFull) == 16, "64-bit clz hand-rolled failed.");
	static_assert(clz(0x000000000FFFFFFFull) == 36, "64-bit clz hand-rolled failed.");
	static_assert(clz(0x0000000000000FFFull) == 52, "64-bit clz hand-rolled failed.");
	static_assert(clz(0x0000000000000000ull) == 64, "64-bit clz hand-rolled failed.");
	static_assert(clz(0x0000000000000003ull) == 62, "64-bit clz hand-rolled failed.");
	static_assert(clz(0x0000013333000003ull) == 23, "64-bit clz hand-rolled failed.");

	/*************clz******************************/

	/*************start detail for estimation***********/
	constexpr double make_alpha(size_t m) {
		switch(m) {
			case 16: return .673;
			case 32: return .697;
			case 64: return .709;
			default: return 0.7213 / (1 + 1.079/m);
		}
	}

	/*************end detail for estimation************/

	static constexpr const char *JESTIM_STRINGS []
	{
		"ORIGINAL", "ERTL_IMPROVED", "ERTL_MLE", "ERTL_JOINT_MLE"
	};
	//	static const char *EST_STRS [] {
	//		"original",
	//			"ertl_improved",
	//			"ertl_mle",
	//			"ertl_joint_mle"
	//	};

	enum EstimationMethod: uint8_t {
		ORIGINAL       = 0,
		ERTL_IMPROVED  = 1, // Improved but biased method
		ERTL_MLE       = 2  // element-wise max, followed by MLE
	};
	enum JointEstimationMethod: uint8_t {
		ERTL_JOINT_MLE = 3  // Ertl special version
	};

	class HyperLogLog{
		protected:
			std::vector<uint8_t> core_;//sketchInfo; TODO: allocator
			mutable double value_; //cardinality
			uint32_t np_; // 10-20
			mutable uint8_t is_calculated_;
			EstimationMethod                        estim_;
			JointEstimationMethod                  jestim_;
			//HashStruct                                 hf_;

		public:
			//HyperLogLog(int np):core_(1uL<<np,0),np_(np),estim_(2),jestim_(3){};
			HyperLogLog(int np):core_(1uL<<np,0),np_(np),is_calculated_(0),estim_(EstimationMethod::ERTL_MLE),jestim_(JointEstimationMethod::ERTL_JOINT_MLE) {};
			~HyperLogLog(){};
			//void update(const char* seq);
			//void showSketch();
			double distance(const HyperLogLog &h2) const {return jaccard_index(h2);}
			//double jaccard_index(HyperLogLog &h2); 
			//double jaccard_index(const HyperLogLog &h2) const; 

			uint32_t p() const {return np_;}
			uint32_t q() const {return (sizeof(uint64_t) * CHAR_BIT) - np_;}
			uint64_t m() const {return static_cast<uint64_t>(1) << np_;}
			size_t size() const {return size_t(m());}
			bool get_is_ready() const {return is_calculated_;}
			const auto &core()    const {return core_;}
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
			//double union_size(const HyperLogLog &other) const;


			//private:
			//void add(uint64_t hashval);
			//void addh(const std::string &element);
			double alpha()          const {return make_alpha(m());}
			static constexpr double LARGE_RANGE_CORRECTION_THRESHOLD = (1ull << 32) / 30.;
			static constexpr double TWO_POW_32 = 1ull << 32;
			static double small_range_correction_threshold(uint64_t m) {return 2.5 * m;}


			// Call sum to recalculate if you have changed contents.
			void sum() const {
				//fprintf(stdout,"run in sum");
				const auto counts(sum_counts(core_)); // std::array<uint32_t, 64>  // K->C
#if DEBUG
				fprintf(stdout,"[W:%s:%d] sum() Counts: ",__PRETTY_FUNCTION__, __LINE__);
				fprintf(stdout,"Counts: ");
				for(int i=0; i<64; i++)
					fprintf(stdout,"%d, ", counts[i]);
				fprintf(stdout,"\n");
#endif
				value_ = calculate_estimate(counts, estim_, m(), np_, alpha(), 1e-2);
				is_calculated_ = 1;
			}
			void csum() const { if(!is_calculated_) sum(); }

			// Returns cardinality estimate. Sums if not calculated yet.
			double creport() const {
				csum();
				return value_;
			}
			double report() noexcept {
				csum();
				return creport();
			}
			//double calculate_estimate(const std::array<uint32_t,64> &counts, EstimationMethod estim, uint64_t m, uint32_t p, double alpha, double relerr) const; 

			//std::array<uint32_t,64> sum_counts(const std::vector<uint8_t> &sketchInfo) const;
			//template<typename T>
			//	void compTwoSketch(const std::vector<uint8_t> &sketch1, const std::vector<uint8_t> &sketch2, T &c1, T &c2, T &cu, T &cg1, T &cg2, T &ceq) const;

			//template<typename FloatType>
			//	static constexpr FloatType gen_sigma(FloatType x); 
			//template<typename FloatType>
			//	static constexpr FloatType gen_tau(FloatType x);
			//template<typename T>
			//	double ertl_ml_estimate(const T& c, unsigned p, unsigned q, double relerr=1e-2) const; 
			//template<typename HllType>
			//	std::array<double, 3> ertl_joint(const HllType &h1, const HllType &h2) const; 

			//	}; // class HyperLogLog


			std::array<uint32_t,64> sum_counts(const std::vector<uint8_t> &sketchInfo) const {
				std::array<uint32_t,64> sum_count {0};//default 64
				for(uint64_t i=0; i<sketchInfo.size(); ++i){
					sum_count[sketchInfo[i]]++;
				}
				return sum_count;
			}

			template<typename T>
				void compTwoSketch(const std::vector<uint8_t> &sketch1, const std::vector<uint8_t> &sketch2, T &c1, T &c2, T &cu, T &cg1, T &cg2, T &ceq) const {
					assert(sketch1.size()==sketch2.size());//
					std::array<uint32_t, 64> c1l{0}, c2l{0}, c1g{0}, c2g{0};
					for(uint64_t i=0; i<sketch1.size(); ++i) {
						//TODO: SIMD
						if(sketch1[i]<sketch2[i]){
							c1l[sketch1[i]]++;
							c2g[sketch2[i]]++;
						} else if(sketch1[i]>sketch2[i]){
							c1g[sketch1[i]]++;
							c2l[sketch2[i]]++;
						} else{
							ceq[sketch1[i]]++;
						}
					}
					for(int i=0; i<64; ++i) {
						c1[i] = c1l[i] + ceq[i] + c1g[i];
						c2[i] = c2l[i] + ceq[i] + c2g[i];
						cu[i] = c1g[i] + ceq[i] + c2g[i];
						cg1[i] = c1g[i];
						cg2[i] = c2g[i];
					}

				}

			// TODO: using SIMD to accelerate
			void update(const char* seq) {
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
					//char kmer[KMERLEN+1];
					char kmer_fwd[KMERLEN+1];
					char kmer_rev[KMERLEN+1];
					memcpy(kmer_fwd, seq+i, KMERLEN);
					memcpy(kmer_rev, seqRev+LENGTH-i-KMERLEN, KMERLEN);
					kmer_fwd[KMERLEN] = '\0';
					kmer_rev[KMERLEN] = '\0';
					if(memcmp(kmer_fwd, kmer_rev, KMERLEN) <= 0) {
						//fprintf(stderr, "kmer_fwd = %s \n", kmer_fwd);
						addh(kmer_fwd);
					} else {
						//fprintf(stderr, "kmer_rev = %s \n", kmer_rev);
						addh(kmer_rev);
					}

				}
			}

			HyperLogLog merge(const HyperLogLog &other) const {
				if(other.p() != p())
					throw std::runtime_error(std::string("p (") + std::to_string(p()) + " != other.p (" + std::to_string(other.p()));
				HyperLogLog ret(*this);
				//ret += other;
				//ret.core_ = max(core_, other.core_);
				for(uint64_t i=0; i<m(); ++i){
					ret.core_[i] = std::max(core_[i],other.core_[i]); 
				}
				return ret;
			}

			//TODO: int hash
			//HyperLogLog::inline void addh(uint64_t element) {
			//	element = hash(element); //TODO: which hf_
			//	add(element);
			//}
			inline void addh(const std::string &element) {
				add(std::hash<std::string>{}(element));
			}

			//TODO: different hash function
			//hash() {
			//}

			//TODO: clz is a function in integral.h
			inline void add(uint64_t hashval) {
				const uint32_t index(hashval >> q());
				const uint8_t lzt(clz(((hashval << 1)|1) << (np_ - 1)) + 1);
				core_[index] = std::max(core_[index], lzt);
#if LZ_COUNTER
				++clz_counts_[clz(((hashval << 1)|1) << (np_ - 1)) + 1];
#endif
			}
			//TODO:Added by liumy to show sketch for testing. 
			void showSketch(){
				int vecSize = core_.size();
				for(int i=0; i<vecSize; ++i)
					fprintf( stdout,"%u, ",  core_[i] );
				fprintf(stdout,"\n");
			}

			double union_size(const HyperLogLog &other) const {
				if(jestim_ != JointEstimationMethod::ERTL_JOINT_MLE) {
					assert(m() == other.m()|| !std::fprintf(stderr, "sizes don't match! Size1: %zu. Size2: %zu\n", m(), other.m()));
					std::array<uint32_t,64> counts{0};
					std::vector<uint8_t> unionCore(m(),0);
					for(uint64_t i=0; i<m(); ++i){
						unionCore[i] = std::max(core_[i],other.core_[i]); 
					}
					counts = sum_counts(unionCore);
					return calculate_estimate(counts, get_estim(), m(), p(), alpha(), 1e-2);
				}
				//std::fprintf(stderr, "jestim is ERTL_JOINT_MLE: %s\n", JESTIM_STRINGS[jestim_]);
				const auto full_counts = ertl_joint(*this, other);
				return full_counts[0] + full_counts[1] + full_counts[2];
			}


			double jaccard_index(const HyperLogLog &h2) const {
				if(jestim_ == JointEstimationMethod::ERTL_JOINT_MLE) {
					auto full_cmps = ertl_joint(*this, h2);
					const auto ret = full_cmps[2] / (full_cmps[0] + full_cmps[1] + full_cmps[2]);
					return ret;
				}
				const double us = union_size(h2);
				const double ret = (creport() + h2.creport() - us) / us;
				return std::max(0., ret);
			}

			//template<typename CountArrType>
			double calculate_estimate(const std::array<uint32_t,64> &counts,
					EstimationMethod estim, uint64_t m, uint32_t p, double alpha, double relerr) const {
				assert(estim <= 3);
#if DEBUG
				fprintf(stdout,"[W:%s:%d] Counts: ",__PRETTY_FUNCTION__, __LINE__);
				fprintf(stdout,"Counts: ");
				for(int i=0; i<64; i++)
					fprintf(stdout,"%d, ", counts[i]);
				fprintf(stdout,"\n");
#endif
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
										   //#if !NDEBUG
										   //										   std::fprintf(stderr, "[W:%s:%d] Small value correction. Original estimate %lf. New estimate %lf.\n",
										   //												   __PRETTY_FUNCTION__, __LINE__, value, m * std::log(static_cast<double>(m) / counts[0]));
										   //#endif
										   value = m * std::log(static_cast<double>(m) / counts[0]);
									   }
								   } else if(value > LARGE_RANGE_CORRECTION_THRESHOLD) {
									   // Reuse sum variable to hold correction.
									   // I do think I've seen worse accuracy with the large range correction, but I would need to rerun experiments to be sure.
									   sum = -std::pow(2.0L, 32) * std::log1p(-std::ldexp(value, -32));
									   if(!std::isnan(sum)) value = sum;
									   //#if !NDEBUG
									   //									   else std::fprintf(stderr, "[W:%s:%d] Large range correction returned nan. Defaulting to regular calculation.\n", __PRETTY_FUNCTION__, __LINE__);
									   //#endif
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

			// Based off https://github.com/oertl/hyperloglog-sketch-estimation-paper/blob/master/c%2B%2B/cardinality_estimation.hpp
			template<typename FloatType>
				static constexpr FloatType gen_sigma (FloatType x) {
					if(x == 1.) return std::numeric_limits<FloatType>::infinity();
					FloatType z(x);
					for(FloatType zp(0.), y(1.); z != zp;) {
						x *= x; zp = z; z += x * y; y += y;
						if(std::isnan(z)) {
							std::fprintf(stderr, "[W:%s:%d] Reached nan. Returning the last usable number.\n", __PRETTY_FUNCTION__, __LINE__);
							return zp;
						}
					}
					return z;
				}

			template<typename FloatType>
				static constexpr FloatType gen_tau (FloatType x) {
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
			template<typename T>
				double ertl_ml_estimate(const T& c, unsigned p, unsigned q, double relerr) const {
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

#if DEBUG
					fprintf(stdout,"[W:%s:%d] Counts: ",__PRETTY_FUNCTION__, __LINE__);
					for(int i=0; i<64; i++)
						fprintf(stdout,"%d, ", c[i]);
					fprintf(stdout,"\n");
#endif

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
					//#if DEBUG
					//					fprintf(stderr,"x = %lf, m = %lf \n", x, m);
					//#endif
					return x*m;
				}
			template<typename HllType>
				std::array<double, 3> ertl_joint(const HllType &h1, const HllType &h2) const {
					assert(h1.m() == h2.m() || !std::fprintf(stderr, "sizes don't match! Size1: %zu. Size2: %zu\n", h1.size(), h2.size()));
					std::array<double, 3> ret;
					if(h1.get_jestim() != JointEstimationMethod::ERTL_JOINT_MLE) {
						// intersection & union
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
					//TODO: K->C5
					//joint_unroller ju;
					//ju.sum_arrays(h1.core(), h2.core(), c1, c2, cu, cg1, cg2, ceq);
					compTwoSketch(h1.core(), h2.core(), c1, c2, cu, cg1, cg2, ceq);
					const double cAX = h1.get_is_ready() ? h1.creport() : ertl_ml_estimate(c1, h1.p(), h1.q(), 1e-2);
					const double cBX = h2.get_is_ready() ? h2.creport() : ertl_ml_estimate(c2, h2.p(), h2.q(), 1e-2);
					const double cABX = ertl_ml_estimate(cu, h1.p(), h1.q(), 1e-2);
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
					double cAXBhalf = ertl_ml_estimate(countsAXBhalf, p, q - 1, 1e-2);
					double cBXAhalf = ertl_ml_estimate(countsBXAhalf, p, q - 1, 1e-2);
					//std::fprintf(stderr, "Made halves: %lf, %lf\n", cAXBhalf, cBXAhalf);
					ret[0] = cABX - cBX;
					ret[1] = cABX - cAX;
					double cX1 = (1.5 * cBX + 1.5*cAX - cBXAhalf - cAXBhalf);
					double cX2 = 2.*(cBXAhalf + cAXBhalf) - 3.*cABX;
					ret[2] = std::max(0., 0.5 * (cX1 + cX2));
					return ret;
				}




		}; // class HyperLogLog

			} //namespace Sketch

#endif // #ifndef HYPERLOGLOG_H_
