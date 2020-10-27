#ifndef Sketch_h
#define Sketch_h

#include <map>
#include <vector>
#include <string>
//#include <string.h>
#include <stdint.h>

#include "MinHash.h"
#include "histoSketch.h"
#include "HyperLogLog.h"

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif  

/// \brief Sketch namespace
namespace Sketch{

	typedef uint64_t hash_t;

	struct WMHParameters
	{
		int kmerSize;
		int sketchSize;
		int windowSize;
		double * r; 
		double * c; 
		double * b; 
	};

	struct Reference
	{
		// no sequence for now
		std::string name;
		std::string comment;
		uint64_t length;
		HashList hashesSorted;
		std::vector<uint32_t> counts;
	};

/// Sketching seqeunces using minhash method
	class MinHash
	{

		public:
			/// minhash init with parameters
			MinHash(int k = 21, int size = 1000, uint32_t seed = 42, bool rc = true):
				kmerSize(k), sketchSize(size), seed(seed), use64(true), noncanonical(!rc)
			{
				minHashHeap = new MinHashHeap(use64, sketchSize);

				this->kmerSpace = pow(alphabetSize, kmerSize);

				this->totalLength = 0;

				this->needToList = true;
			}

			/// minhash is updatable with multiple sequences
			void update(char * seq);

			/// merge two minhashes
			void merge(MinHash& msh);

			/// return the jaccard index
			double jaccard(MinHash * msh);			

			/// return mutation distance defined in Mash instead of jaccard distance
			double mdistance(MinHash * msh);

			//
			//void getMinHash(); //return hash values to a vector?

			//print hash values for debug
			void printMinHashes();

			/// return totalSeqence length, including multiple updates
			uint64_t getTotalLength(){return totalLength;}

			///Estimate the cardinality count
			int count(){ return minHashHeap->estimateSetSize();}
			// parameters

			/// return kmerSize
			int getKmerSize() { return kmerSize; }

			// return alphabet size
			//uint32_t getAlphabetSize() { return alphabetSize; }

			// return whether to preserve case
			//bool isPreserveCase() { return preserveCase; }

			// return whether to use 64bit hash
			//bool isUse64() { return use64; }

			/// return hash seed
			uint32_t getSeed() {return seed; }

			/// return sketch size
			uint32_t getMaxSketchSize() {return sketchSize; }

			/// return whether to use reverse complement
			bool isReverseComplement() { return !noncanonical; }

			/// test whether this minhash is empty
			bool isEmpty() { 
				if(this->needToList){
					this->heapToList();
					this->needToList = false;
				}

				if(this->reference.hashesSorted.size() <= 0)
					return true;
				else
					return false;
				
			}

			/// get sketch size, it should be less than max sketch size
			int getSketchSize() {
				if(this->needToList)
				{
					this->heapToList();
					this->needToList = false;
				}
				
				return this->reference.hashesSorted.size();
			}

		private:
			bool needToList = true;
			double kmerSpace;
			MinHashHeap * minHashHeap;
			Reference reference;
			uint64_t totalLength = 0;

			double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize);
			void heapToList();

			//parameters
			int kmerSize;
			uint32_t alphabetSize; //nuc sequences
			uint32_t seed;
			uint64_t sketchSize; //minHashesPerWindow
			bool noncanonical;
			bool use64; //always true (remove support for hash32)

			//FIXME: perserveCase is not included in constructor
			bool preserveCase = false;
	};

	class WMinHash{
		public:
			//WMinHash(int k = 21, int size = 50, int windowSize = 20, double paraDWeight = 0.0):
			WMinHash(WMHParameters parameters, double paraDWeight = 0.0):
				kmerSize(parameters.kmerSize), histoSketchSize(parameters.sketchSize), minimizerWindowSize(parameters.windowSize), paraDecayWeight(paraDWeight)
			{	
				//numBins = pow(kmerSize, alphabetSize); //need to be confirmed
				//histoDimension = pow(kmerSize, alphabetSize); //need to be confirmed
				numBins = pow(kmerSize, 4); //consult from HULK
				histoDimension = numBins; //need to be confirmed
			
				binsArr = (double *)malloc(numBins * sizeof(double));
				memset(binsArr, 0, numBins*sizeof(double));//init binsArr
			
				int g = ceil(2 / EPSILON);
				int d = ceil(log(1 - DELTA) / log(0.5));
				countMinSketch = (double *)malloc(d * g *sizeof(double));
				memset(countMinSketch, 0, d*g*sizeof(double));//init countMinSketch 

				r = parameters.r;
				c = parameters.c;
				b = parameters.b;
				
			//	//the r, c, b and getCWS need to be outClass
			//	r = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
			//	c = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
			//	b = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
			//	getCWS(r, c, b, histoSketchSize, histoDimension);
			
				////for getCWS test debug
			//	FILE * fp;
			//	fp = fopen("cws.txt", "w");
			//	for(int i = 0; i < histoSketchSize*histoDimension; i++){
			//		fprintf(fp, "%ld\t%lf\t%lf\t%lf\n", i, r[i], c[i], b[i]);
			//	}
			
				histoSketches = (uint32_t *) malloc (histoSketchSize * sizeof(uint32_t));
				histoWeight = (double *) malloc (histoSketchSize * sizeof(double));
				
				//add the applyConceptDrift and decayWeight.
				if(paraDecayWeight < 0.0 || paraDecayWeight > 1.0){
					cerr << "the paraDecayWeight must between 0.0 and 1.0 " << endl;
					exit(1);
				}
				else{
					applyConceptDrift = true;
				}
				if(paraDecayWeight == 1.0){
					applyConceptDrift = false;
				}
				decayWeight = 1.0;
				if(applyConceptDrift){
					decayWeight = exp(-paraDecayWeight);
				}
			
				needToCompute = true;
			
			}
			

			~WMinHash();

			/// weightedMinHash is updatable with multiple sequences
			void update(char * seq);

			/// return the jaccard index
			double wJaccard(WMinHash * wmh);

			/// return the distance which is negative correlation with jaccard index
			double distance(WMinHash * wmh);

			/// print hash value fo debug
			void getWMinHash();

			//void setKmerSize(int kmerSizeNew) { kmerSize = kmerSizeNew; }
			//void setAlphabetSize(int alphabetSizeNew) { alphabetSize = alphabetSizeNew; }
			//void setNumBins(int numBinsNew) { numBins = numBinsNew; }
			//void setMinimizerWindowSize(int minimizerWindowSizeNew) { minimizerWindowSize = minimizerWindowSizeNew; }
			//void setHistoSketchSize(int histoSketchSizeNew); //{ histoSketchSize = histoSketchSizeNew; }
			//void setHistoDimension(int histoDimensionNew); //{ histoDimension = histoDimensionNew; }
			//void setParaDecayWeight(double paraDecayWeightNew) { paraDecayWeight = paraDecayWeightNew; }
			//void setApplyConceptDrift(bool applyConceptDriftNew) { applyConceptDrift = applyConceptDriftNew; }

			/// return kmerSize
			int getKmerSize() { return kmerSize; }

			/// return the minimizer window size
			int getMinimizerWindowSize() { return minimizerWindowSize; }

			//int getAlphabetSize() { return alphabetSize; }
			//int getNumBins() { return numBins; }
			//int getMinimizerWindowSize() { return minimizerWindowSize; }
			//int getHistoSketchSize() { return histoSketchSize; }
			//int getHistoDimension() { return histoDimension; }
			//double getParaDecayWeight() { return paraDecayWeight; }
			//bool isApplyComceptDrift() { return applyConceptDrift; }


			

		private:
			WMHParameters parameters;
			int kmerSize;
			int alphabetSize = 4;
			int numBins;
			int minimizerWindowSize;
			int histoSketchSize;
			int histoDimension;
			double paraDecayWeight;
			bool applyConceptDrift;
			double decayWeight;

			bool needToCompute = true;
			double * binsArr;
			double * countMinSketch; 
			std::vector<uint64_t> sketches;
			std::vector<Bin> kmerSpectrums;
			uint32_t * histoSketches;
			double * histoWeight;

			void computeHistoSketch();

			double * r;
			double * c;
			double * b;



	};

	//OMinHash
	struct OSketch {
		//std::string       name;
		int               k, l, m;
		std::vector<char> data;
		std::vector<char> rcdata;

		bool operator==(const OSketch& rhs) const {
			return k == rhs.k && l == rhs.l && m == rhs.m && data == rhs.data && rcdata == rhs.rcdata;
		}

	};

	/// Sketching and compare sequences or strings using Order MinHash algorithm.
	class OrderMinHash{

		public:
			/// OrderMinHash constructor
			OrderMinHash(){};
			/// OrderMinHash constructor for sketching sequences using default parameters
			OrderMinHash(char * seqNew);
			~OrderMinHash() {if (rcseq != NULL) delete rcseq;};

			/// return sketch result in `OSketch` type
			OSketch getSektch(){ return sk;}

			/** \rst
			  Build a `OrderMinHash` sketch.
			  `seqNew` is NULL pointer in default.
			  If seqNew is NULL pointer, buildSketch() will rebuild sketh using old data.
			  This is useful when chaning parameters and build a new sketch.
			 \endrst
			*/
			void buildSketch(char * seqNew);

			/** 
			   Return similarity between two `OrderMinHash` sketches. In `OrderMinHash` class, there is no jaccard function provied. Because `OrderMinHash` is a proxy of edit distance instead of jaccard index.
			*/
			double similarity(OrderMinHash & omh2);

			/// Return distance between two `OrderMinHash` sketches
			double distance(OrderMinHash & omh2)
			{
				return (double)1.0 - similarity(omh2);
			}

			/// Set parameter `kmerSize`: default 21.
			void setK(int k){ m_k = k; }

			/// Set parameter `l`: default 2 (normally 2 - 5).
			void setL(int l){ m_l = l; }

			/// Set parameter `m`: default 500.
			void setM(int m){ m_m = m; }

			/// Set seed value for random generator: default 32.
			void setSeed(uint64_t seedNew) { mtSeed = seedNew; }

			/** 
			  Choose whether to deal with reverse complement sequences: default false.
			  Reverse complement is normally used in biological sequences such as DNA or protein sequences.
			*/
			void setReverseComplement(bool isRC){rc = isRC;}

			/// Return parameter `kmerSize`.
			int getK(){return m_k;}	

			/// Return parameter `l`.
			int getL(){return m_l;}	

			/// Return parameter `m`.
			int getM(){return m_m;}		

			/// Return random generator seed value.
			uint64_t getSeed() { return mtSeed; }

			/// Test whether to deal with reverse complement kmers.
			bool isReverseComplement(){return rc;}

		private:

			char * seq = NULL;
			char * rcseq = NULL;
			//Parameters parameters;
			int m_k = 21, m_l = 2, m_m = 500;
			//reverse complement
			bool rc = false;
			OSketch sk;
			uint64_t mtSeed = 32; //default value

			void sketch();

			inline void compute_sketch(char * ptr, const char * seq);

			double compare_sketches(const OSketch& sk1, const OSketch& sk2, 
											  ssize_t m = -1, bool circular = false);
			double compare_sketch_pair(const char* p1, const char* p2,
									   unsigned m, unsigned k, unsigned l, bool circular);

	};
	
	class HyperLogLog{
		
		public:
			HyperLogLog(int np):core_(1uL<<np,0),np_(np),is_calculated_(0),estim_(EstimationMethod::ERTL_MLE),jestim_(JointEstimationMethod::ERTL_JOINT_MLE) {};
			~HyperLogLog(){};
			void update(const char* seq);
			HyperLogLog merge(const HyperLogLog &other) const;
			void printSketch();
			double distance(const HyperLogLog &h2) const {return jaccard_index(h2);}
			double jaccard_index(HyperLogLog &h2); 
			double jaccard_index(const HyperLogLog &h2) const; 

		protected:
			std::vector<uint8_t> core_;//sketchInfo; 
			mutable double value_; //cardinality
			uint32_t np_; // 10-20
			mutable uint8_t is_calculated_;
			EstimationMethod                        estim_;
			JointEstimationMethod                  jestim_;
			//HashStruct                                 hf_;
		
		private:
			uint32_t p() const {return np_;}//verification
			uint32_t q() const {return (sizeof(uint64_t) * CHAR_BIT) - np_;}
			uint64_t m() const {return static_cast<uint64_t>(1) << np_;}
			size_t size() const {return size_t(m());}
			bool get_is_ready() const {return is_calculated_;}
			const auto &core()    const {return core_;}
			EstimationMethod get_estim()       const {return  estim_;}
			JointEstimationMethod get_jestim() const {return jestim_;}
			void set_estim(EstimationMethod val) { estim_ = std::max(val, ERTL_MLE);}
			void set_jestim(JointEstimationMethod val) { jestim_ = val;}
			void set_jestim(uint16_t val) {set_jestim(static_cast<JointEstimationMethod>(val));}
			void set_estim(uint16_t val)  {estim_  = static_cast<EstimationMethod>(val);}

			// Returns cardinality estimate. Sums if not calculated yet.
			double creport() const {
				csum();
				return value_;
			}
			double report() noexcept {
				csum();
				return creport();
			}


			//private:
			void add(uint64_t hashval);
			void addh(const std::string &element);
			double alpha()          const {return make_alpha(m());}
			static double small_range_correction_threshold(uint64_t m) {return 2.5 * m;}
			double union_size(const HyperLogLog &other) const;
			// Call sum to recalculate if you have changed contents.
			void csum() const { if(!is_calculated_) sum(); }
			void sum() const {
				const auto counts(sum_counts(core_)); // std::array<uint32_t, 64>  // K->C
				value_ = calculate_estimate(counts, estim_, m(), np_, alpha(), 1e-2);
				is_calculated_ = 1;
			}
			std::array<uint32_t,64> sum_counts(const std::vector<uint8_t> &sketchInfo) const;
			double calculate_estimate(const std::array<uint32_t,64> &counts, EstimationMethod estim, uint64_t m, uint32_t p, double alpha, double relerr) const; 
			template<typename T>
				void compTwoSketch(const std::vector<uint8_t> &sketch1, const std::vector<uint8_t> &sketch2, T &c1, T &c2, T &cu, T &cg1, T &cg2, T &ceq) const;
			template<typename T>
				double ertl_ml_estimate(const T& c, unsigned p, unsigned q, double relerr=1e-2) const; 
			template<typename HllType>
				std::array<double, 3> ertl_joint(const HllType &h1, const HllType &h2) const; 



	};


}//namespace sketch

#endif //Sketch_h
