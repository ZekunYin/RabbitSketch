#ifndef Sketch_h
#define Sketch_h

#include <map>
#include <vector>
#include <string>
//#include <string.h>
#include <stdint.h>

#include "MinHash.h"
#include "histoSketch.h"
#include "hll/hyperloglog.h"


/// \brief Sketch namespace
namespace Sketch{

	typedef uint64_t hash_t;

	class Parameters
	{
	public:
		Parameters()
			:
				//MinHash
				kmerSize(21),
				alphabetSize(4),
				preserveCase(false),
				use64(true),
				seed(42),
				minHashesPerWindow(1000),
				noncanonical(false),

				//for WMinHash
				numBins(194481),
				minimizerWindowSize(9),
				histoSketch_sketchSize(50),
				histoSketch_dimension(194481),
				//applyConceptDrift(false),
				paraDecayWeight(0.0),

				//for OMinHash
				l(2),
				m(500),
				rc(false)
		{
			memset(alphabet, 0, 256);
			getCWS(r, c, b, 50, 194481);
			
		}

		Parameters(const Parameters & other)
			:
				//MinHash	
				kmerSize(other.kmerSize),
				alphabetSize(other.alphabetSize),
				preserveCase(other.preserveCase),
				use64(other.use64),
				seed(other.seed),
				minHashesPerWindow(other.minHashesPerWindow),
				noncanonical(other.noncanonical),

				//for WMinHash
				numBins(other.numBins),
				minimizerWindowSize(other.minimizerWindowSize),
				histoSketch_sketchSize(other.histoSketch_sketchSize),
				histoSketch_dimension(other.histoSketch_dimension),
				//applyConceptDrift(other.applyConceptDrift),
				paraDecayWeight(other.paraDecayWeight),
				r(other.r),
				c(other.c),
				b(other.b),

				//for OMinHash
				l(other.l),
				m(other.m),
				rc(other.rc)

		{
			memcpy(alphabet, other.alphabet, 256);
		}
//		int get_kmerSize(){ return kmerSize; }
//		uint32_t get_alphabetSize(){ return alphabetSize;}
//		bool get_preserveCase(){ return preserveCase;}
//		bool get_use64(){ return use64;}
		int get_numBins(){return numBins;}
		int get_minimizerWindowSize(){return minimizerWindowSize;}
		int get_histoSketch_sketchSize(){return histoSketch_sketchSize;}
		int get_histoSketch_dimension(){return histoSketch_dimension;}
		double get_paraDecayWeight(){return paraDecayWeight;}
		double * getR(){
			//printf("the point in getR of r is: %p\n",r);
			return r;
		}
		double * getC(){
			//printf("the point in getC of c is: %p\n",c);
			return c;
		}
		double * getB(){
			//printf("the point in getB of b is: %p\n",b);
			return b;
		}
		

		void set_numBins(int input){numBins = input;}
		void set_minimizerWindowSize(int input){minimizerWindowSize = input;}
		void set_histoSketch_sketchSize(int input){
			histoSketch_sketchSize = input;
			r = (double *)malloc(histoSketch_sketchSize * histoSketch_dimension * sizeof(double));
			c = (double *)malloc(histoSketch_sketchSize * histoSketch_dimension * sizeof(double));
			b = (double *)malloc(histoSketch_sketchSize * histoSketch_dimension * sizeof(double));
			getCWS(r, c, b, histoSketch_sketchSize, histoSketch_dimension);
			printf("the point in set histoSketch_sketchSize of r is: %p\n",r);
			printf("the point in set histoSketch_sketchSize of c is: %p\n",c);
			printf("the point in set histoSketch_sketchSize of b is: %p\n",b);
		}
		void set_histoSketch_dimension(int input){
			histoSketch_dimension = input;
			r = (double *)malloc(histoSketch_sketchSize * histoSketch_dimension * sizeof(double));
			c = (double *)malloc(histoSketch_sketchSize * histoSketch_dimension * sizeof(double));
			b = (double *)malloc(histoSketch_sketchSize * histoSketch_dimension * sizeof(double));
			getCWS(r, c, b, histoSketch_sketchSize, histoSketch_dimension);
			printf("the point in set histoSketch_dimension of r is: %p\n",r);
			printf("the point in set histoSketch_dimension of c is: %p\n",c);
			printf("the point in set histoSketch_dimension of b is: %p\n",b);
		}
		void set_paraDecayWeight(double input){paraDecayWeight = input;}

		



//	private:
		//MinHash
		int kmerSize;
		bool alphabet[256];
		uint32_t alphabetSize;
		bool preserveCase;
		bool use64;
		uint32_t seed;
		uint64_t minHashesPerWindow; //SketchSize
		bool noncanonical;

		//parameters for order minhash
		int l;
		int m;
		bool rc;
		//kmerSize for k
	

	private:
		//parameters for weight minHash@xxm
		int numBins;
		int minimizerWindowSize;
		int histoSketch_sketchSize;
		int histoSketch_dimension;
		//bool applyConceptDrift;
		double paraDecayWeight;

		//double * r;
		//double * c;
		//double * b;

		double * r = (double *)malloc(50 * 194481 * sizeof(double));
		double * c = (double *)malloc(50 * 194481 * sizeof(double));
		double * b = (double *)malloc(50 * 194481 * sizeof(double));



	};


	//typedef robin_hood::unordered_set<hash_t> Hash_set;

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
			MinHash(Parameters parametersNew);
			/// minhash is updatable with multiple sequences
			void update(char * seq);
			/// merge two minhashes
			void merge(MinHash& msh);
			/// return the jaccard index
			double jaccard(MinHash * msh);			
			/// return distance defined in Mash and RabbitMash
			double dist(MinHash * msh);
			//HashList & getHashList(); //TODO: return to vector instead of HashList
			void getMinHash(); //return hash values to a vector?
			void printMinHashes();
			uint64_t getTotalLength(){return totalLength;}//return totalSeqence length.

		private:
			bool needToList = true;
			double kmerSpace;
			Parameters parameters;
			MinHashHeap * minHashHeap;
			Reference reference;
			uint64_t totalLength = 0;
			double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize);
			void heapToList();

	};

	class WMinHash{
		public:
			WMinHash(Parameters parametersNew);
			~WMinHash();

			void update(char * seq);

			double wJaccard(WMinHash * wmh);
			double distance(WMinHash * wmh);
			
			void getWMinHash();

		private:
			Parameters parameters;
			bool needToCompute = true;
			double * binsArr;
			double * countMinSketch; 
			std::vector<uint64_t> sketches;
			std::vector<Bin> kmerSpectrums;

			void computeHistoSketch();

			double * r;
			double * c;
			double * b;
			uint32_t * histoSketch_sketch;
			double * histoSketch_sketchWeight;
			bool applyConceptDrift;
			double decayWeight;



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

	class OMinHash{

		public:
			//OMinHash(Parameters parametersNew);
			OMinHash(Parameters parametersNew, char * seqNew);
			~OMinHash() {if (rcseq != NULL) delete rcseq;};
			//update -- not supported yet!

			void sketch();
			OSketch getSektch(){ return sk;}

			double similarity(OMinHash & omh2);

			double distance(OMinHash & omh2){
				return (double)1.0 - similarity(omh2);
			}

		protected:

		private:
			char * seq = NULL;
			char * rcseq = NULL;
			Parameters parameters;
			int m_k, m_l, m_m;
			bool rc;//reverse complement
			OSketch sk;

			inline void compute_sketch(char * ptr, const char * seq);

			double compare_sketches(const OSketch& sk1, const OSketch& sk2, 
											  ssize_t m = -1, bool circular = false);
			double compare_sketch_pair(const char* p1, const char* p2,
									   unsigned m, unsigned k, unsigned l, bool circular);

	};

}//namespace sketch

#endif //Sketch_h
