#include <stdio.h>
#include <iostream>
#include <string.h>
#include <list>
#include <cstring>
#include "Sketch.h"
#include "MurmurHash3.h"
#include "hash.h"
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <immintrin.h>
#include <stdint.h>
#include "xxhash.hpp"
#include <random>
#include "countMin.h" //for EPSILON and DELTA@xxm
#include "minimizer.h"
#include "cws.h"
#include "histoSketch.h"

//#ifdef USE_BOOST
//    #include <boost/math/distributions/binomial.hpp>
//    using namespace::boost::math;
//#else
//    #include <gsl/gsl_cdf.h>
//#endif


using namespace std;
using namespace Sketch;

__m512i inline min512(__m512i v1, __m512i v2)
{
	__mmask8 msk_gt, msk_lt;
	msk_gt = _mm512_cmpgt_epi64_mask(v1, v2);
	msk_lt = _mm512_cmplt_epi64_mask(v1, v2);

	return (msk_gt < msk_lt) ? v1 : v2;
}
void inline transpose8_epi64(__m512i *row0, __m512i* row1, __m512i* row2,__m512i* row3, __m512i* row4, __m512i * row5,__m512i* row6,__m512i * row7)
{
	__m512i __t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7;
	__m512i __tt0,__tt1,__tt2,__tt3,__tt4,__tt5,__tt6,__tt7;

	__m512i idx1,idx2;
	idx1 = _mm512_set_epi64(0xD,0xC,0x5,0x4,0x9,0x8,0x1,0x0);
	idx2 = _mm512_set_epi64(0xF,0xE,0x7,0x6,0xB,0xA,0x3,0x2);

	__t0 = _mm512_unpacklo_epi64(*row0,*row1);
	__t1 = _mm512_unpackhi_epi64(*row0,*row1);
	__t2 = _mm512_unpacklo_epi64(*row2,*row3);
	__t3 = _mm512_unpackhi_epi64(*row2,*row3);
	__t4 = _mm512_unpacklo_epi64(*row4,*row5);
	__t5 = _mm512_unpackhi_epi64(*row4,*row5);
	__t6 = _mm512_unpacklo_epi64(*row6,*row7);
	__t7 = _mm512_unpackhi_epi64(*row6,*row7);


	__tt0 = _mm512_permutex2var_epi64(__t0,idx1,__t2);
	__tt2 = _mm512_permutex2var_epi64(__t0,idx2,__t2);
	__tt1 = _mm512_permutex2var_epi64(__t1,idx1,__t3);
	__tt3 = _mm512_permutex2var_epi64(__t1,idx2,__t3);
	__tt4 = _mm512_permutex2var_epi64(__t4,idx1,__t6);
	__tt6 = _mm512_permutex2var_epi64(__t4,idx2,__t6);
	__tt5 = _mm512_permutex2var_epi64(__t5,idx1,__t7);
	__tt7 = _mm512_permutex2var_epi64(__t5,idx2,__t7);

	*row0 = _mm512_shuffle_i64x2(__tt0,__tt4,0x44);
	*row1 = _mm512_shuffle_i64x2(__tt1,__tt5,0x44);
	*row2 = _mm512_shuffle_i64x2(__tt2,__tt6,0x44);
	*row3 = _mm512_shuffle_i64x2(__tt3,__tt7,0x44);

	*row4 = _mm512_shuffle_i64x2(__tt0,__tt4,0xEE);
	*row5 = _mm512_shuffle_i64x2(__tt1,__tt5,0xEE);
	*row6 = _mm512_shuffle_i64x2(__tt2,__tt6,0xEE);
	*row7 = _mm512_shuffle_i64x2(__tt3,__tt7,0xEE);
}

//Writen by qzh.
//Sketch::MinHash::MinHash(const Sketch::Parameters & parameters)
MinHash ::MinHash(Parameters parametersNew):parameters(parametersNew)
{
	minHashHeap = new MinHashHeap(parameters.use64, parameters.minHashesPerWindow, parameters.reads ?  parameters.minCov : 1);

	this->kmerSpace = pow(parameters.alphabetSize, parameters.kmerSize);
	//cerr << "kmerSpace init from pow is " << this->kmerSpace << endl;
	this->length = 0;
	this->needToList = true;

}

void MinHash::update(char * seq)
{
	//addMinHashes(minHashHeap, seq, LENGTH, parameters);
	//cout << "seq: " << seq << endl;
	const uint64_t LENGTH = strlen(seq);
	this->length += LENGTH;
	int kmerSize = parameters.kmerSize;
	uint64_t mins = parameters.minHashesPerWindow;
	bool noncanonical = parameters.noncanonical;//False.

	// uppercase TODO: alphabets?
	for ( uint64_t i = 0; i < LENGTH; i++ )
	{
		if ( ! parameters.preserveCase && seq[i] > 96 && seq[i] < 123 )
		{
			seq[i] -= 32;
		}
	}

	char * seqRev;

	if ( ! noncanonical )
	{
		seqRev = new char[LENGTH];
		//reverseComplement(seq, seqRev, length);

		char table[4] = {'T','G','A','C'};
		for ( uint64_t i = 0; i < LENGTH; i++ )
		{
			char base = seq[i];

			base >>= 1;
			base &= 0x03;
			seqRev[LENGTH - i - 1] = table[base];


			//		    switch ( base )
			//		    {
			//		        case 'A': base = 'T'; break;
			//		        case 'C': base = 'G'; break;
			//		        case 'G': base = 'C'; break;
			//		        case 'T': base = 'A'; break;
			//		        default: break;
			//		    }
			//		    seqRev[LENGTH - i - 1] = base;
		}
	}

	//for ( uint64_t i = 0; i < LENGTH - kmerSize + 1; i++ )
	//{
	//	// repeatedly skip kmers with bad characters
	//	bool bad = false;
	//	
	//	//Modified by qzh.To detect the correct of alphabet, but it consumes too much time. So we should optimize this process.
	//	//for ( uint64_t j = i; j < i + kmerSize && i + kmerSize <= LENGTH; j++ )
	//	//{
	//	//	if ( ! parameters.alphabet[seq[j]] )
	//	//	{
	//	//		i = j; // skip to past the bad character
	//	//		bad = true;
	//	//		break;
	//	//	}
	//	//}
	//	//
	//	if ( bad )
	//	{
	//		continue;
	//	}
	//	//	
	//	if ( i + kmerSize > LENGTH )
	//	{
	//		// skipped to end
	//		break;
	//	}
	//        
	//    const char *kmer_fwd = seq + i;
	//    const char *kmer_rev = seqRev + LENGTH - i - kmerSize;
	//    const char * kmer = (noncanonical || memcmp(kmer_fwd, kmer_rev, kmerSize) <= 0) ? kmer_fwd : kmer_rev;
	//    bool filter = false;
	//    
	//    hash_u hash = getHash(kmer, kmerSize, parameters.seed, parameters.use64);
	//    
	//	minHashHeap -> tryInsert(hash);
	//}

	//============================================================================================================================================================
	//cout << "use the intel avx512 " << endl;
	//exit(0);
	//addbyxxm
	//hash_u hash = getHash(kmer, kmerSize, parameters.seed, parameters.use64);
	//  const char * kmer = (noncanonical || memcmp(kmer_fwd, kmer_rev, kmerSize) <= 0) ? kmer_fwd : kmer_rev;
	int pend_k = ((kmerSize - 1) / 16 + 1) * 16;
	int n_kmers = length - kmerSize + 1;
	int n_kmers_body = (n_kmers / 16) * 16;
	//int n_kmers_body = (n_kmers / 32) * 32;
	//int n_kmers_body = (n_kmers / 8) * 8;

	const uint8_t* input8 = (const uint8_t *)seq;
	const uint8_t* input8_rev = (const uint8_t *)seqRev;
	//uint64_t* res = (uint64_t* )_mm_malloc(n_kmers * 2 * sizeof(uint64_t), 64);
	uint64_t res[16 * 2];
	//uint64_t res[32 * 2];
	//uint64_t res[8 * 2];
	uint64_t res2[2];
	uint8_t kmer_buf[kmerSize];
	//	uint64_t result[16];

	__m512i v0, v1;
	__m512i vi[8];
	__m512i vj[8];
	//	__m512i vk[8];
	//	__m512i vl[8];
	__m512i vzero = _mm512_setzero_si512();

	__mmask64 mask_load = 0xffffffffffffffff;
	mask_load >>= (64 - kmerSize);

	for(int i = 0; i < n_kmers_body-1; i+=16){

		vi[0] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 0, input8_rev + length - i - 0 - kmerSize, kmerSize) <= 0 ? input8 + i + 0 : input8_rev + length - i - 0 - kmerSize);
		vi[1] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 1, input8_rev + length - i - 1 - kmerSize, kmerSize) <= 0 ? input8 + i + 1 : input8_rev + length - i - 1 - kmerSize);
		vi[2] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 2, input8_rev + length - i - 2 - kmerSize, kmerSize) <= 0 ? input8 + i + 2 : input8_rev + length - i - 2 - kmerSize);
		vi[3] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 3, input8_rev + length - i - 3 - kmerSize, kmerSize) <= 0 ? input8 + i + 3 : input8_rev + length - i - 3 - kmerSize);
		vi[4] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 4, input8_rev + length - i - 4 - kmerSize, kmerSize) <= 0 ? input8 + i + 4 : input8_rev + length - i - 4 - kmerSize);
		vi[5] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 5, input8_rev + length - i - 5 - kmerSize, kmerSize) <= 0 ? input8 + i + 5 : input8_rev + length - i - 5 - kmerSize);
		vi[6] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 6, input8_rev + length - i - 6 - kmerSize, kmerSize) <= 0 ? input8 + i + 6 : input8_rev + length - i - 6 - kmerSize);
		vi[7] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 7, input8_rev + length - i - 7 - kmerSize, kmerSize) <= 0 ? input8 + i + 7 : input8_rev + length - i - 7 - kmerSize);

		vj[0] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 8, input8_rev + length - i - 8 - kmerSize, kmerSize) <= 0 ? input8 + i + 8 : input8_rev + length - i - 8 - kmerSize);
		vj[1] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 9, input8_rev + length - i - 9 - kmerSize, kmerSize) <= 0 ? input8 + i + 9 : input8_rev + length - i - 9 - kmerSize);
		vj[2] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 10, input8_rev + length - i - 10 - kmerSize, kmerSize) <= 0 ? input8 + i + 10 : input8_rev + length - i - 10 - kmerSize);
		vj[3] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 11, input8_rev + length - i - 11 - kmerSize, kmerSize) <= 0 ? input8 + i + 11 : input8_rev + length - i - 11 - kmerSize);
		vj[4] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 12, input8_rev + length - i - 12 - kmerSize, kmerSize) <= 0 ? input8 + i + 12 : input8_rev + length - i - 12 - kmerSize);
		vj[5] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 13, input8_rev + length - i - 13 - kmerSize, kmerSize) <= 0 ? input8 + i + 13 : input8_rev + length - i - 13 - kmerSize);
		vj[6] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 14, input8_rev + length - i - 14 - kmerSize, kmerSize) <= 0 ? input8 + i + 14 : input8_rev + length - i - 14 - kmerSize);
		vj[7] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 15, input8_rev + length - i - 15 - kmerSize, kmerSize) <= 0 ? input8 + i + 15 : input8_rev + length - i - 15 - kmerSize);



		transpose8_epi64(&vi[0], &vi[1], &vi[2], &vi[3], &vi[4], &vi[5], &vi[6], &vi[7]); 
		transpose8_epi64(&vj[0], &vj[1], &vj[2], &vj[3], &vj[4], &vj[5], &vj[6], &vj[7]); 
		//	transpose8_epi64(&vk[0], &vk[1], &vk[2], &vk[3], &vk[4], &vk[5], &vk[6], &vk[7]); 
		//	transpose8_epi64(&vl[0], &vl[1], &vl[2], &vl[3], &vl[4], &vl[5], &vl[6], &vl[7]); 

		//MurmurHash3_x64_128_avx512_8x16(vi, vj, pend_k, kmerSize, 42, &res[2 * i]);// the seed in Mash is 42; verified by xxm;
		MurmurHash3_x64_128_avx512_8x16(vi, vj, pend_k, kmerSize, 42, res);// the seed in Mash is 42; verified by xxm;
		//MurmurHash3_x64_128_avx512_8x32(vi, vj, vk, vl, pend_k, kmerSize, 42, res);// the seed in Mash is 42; verified by xxm;
		//MurmurHash3_x64_128_avx512_8x8(vi, pend_k, kmerSize, 42, res);// the seed in Mash is 42; verified by xxm;

		hash_u hash;
		for(int j = 0; j < 16; j++){
			hash.hash64 = res[j * 2];
			//cout << hash.hash64 << endl;
			minHashHeap->tryInsert(hash);
		}
	}

	for(int i = n_kmers_body; i < n_kmers; i++){
		bool noRev = (memcmp(input8 + i, input8_rev + length - i - kmerSize, kmerSize) <= 0);
		if(noRev){
			for(int j = 0; j < kmerSize; j++){
				//kmer_buf[j] = input8[i + j];
				kmer_buf[j] = input8[i + j];
			}
		}
		else{
			for(int j = 0; j < kmerSize; j++){
				//kmer_buf[j] = input8[i + j];
				kmer_buf[j] = input8_rev[length - i - kmerSize + j];
			}
		}


		//MurmurHash3_x64_128(kmer_buf, kmerSize, 42, &res[2 * i]);// the getHash just need the lower 64bit of the total 128bit of res[i];
		MurmurHash3_x64_128(kmer_buf, kmerSize, 42, res2);// the getHash just need the lower 64bit of the total 128bit of res[i];
		hash_u hash;
		hash.hash64 = res2[0];
		minHashHeap->tryInsert(hash);

	}
	//============================================================================================================================================================

	if ( ! noncanonical )
	{
		delete [] seqRev;
	}
	//	if(needToList){
	//		heapToList();
	//	}
	needToList = true;
}

void MinHash::heapToList()
{
	HashList & hashlist = reference.hashesSorted;
	hashlist.clear();
	hashlist.setUse64(parameters.use64);
	minHashHeap -> toHashList(hashlist);
	minHashHeap -> toCounts(reference.counts);
	hashlist.sort();

}

void MinHash::getMinHash()
{
	//setMinHashesForReference(reference, minHashHeap);
	//	HashList & hashlist = reference.hashesSorted;
	//	hashlist.clear();
	//	hashlist.setUse64(parameters.use64);
	//	minHashHeap -> toHashList(hashlist);
	//	minHashHeap -> toCounts(reference.counts);
	//	hashlist.sort();
	if(needToList){
		heapToList();
		needToList = false;
	}

	for(int i = 0; i < reference.hashesSorted.size(); i++){
		if(parameters.use64)
			cerr << "hash64 " <<  i << " " << reference.hashesSorted.at(i).hash64 << endl;
		else
			cerr << "hash32 " <<  i << " " << reference.hashesSorted.at(i).hash32 << endl;
	}
	return;
}

hash_u getHash(const char * seq, int length, uint32_t seed, bool use64)
{

#ifdef ARCH_32
	char data[use64 ? 8 : 4];
	MurmurHash3_x86_32(seq, length > 16 ? 16 : length, seed, data);
	if ( use64 )
	{
		MurmurHash3_x86_32(seq + 16, length - 16, seed, data + 4);
	}
#else
	char data[16];
	MurmurHash3_x64_128(seq, length, seed, data);
#endif

	hash_u hash;

	if ( use64 )
	{
		hash.hash64 = *((hash64_t *)data);
	}
	else
	{
		hash.hash32 = *((hash32_t *)data);
	}

	return hash;
}

bool hashLessThan(hash_u hash1, hash_u hash2, bool use64)
{
	if ( use64 )
	{
		return hash1.hash64 < hash2.hash64;
	}
	else
	{
		return hash1.hash32 < hash2.hash32;
	}
}

double MinHash::jaccard(MinHash * msh)
{
	if(needToList){
		heapToList();
		needToList = false;
	}
	if(msh->needToList){
		cout << "msh2 need to list addbyxxm " << endl;
		msh->heapToList();
		msh->needToList = false;
	}

	//transport the previous parameter @xxm
	uint64_t sketchSize = this->parameters.minHashesPerWindow;
	int kmerSize = this->parameters.kmerSize;
	//int kmerSpace = this->kmerSpace;

	uint64_t i = 0;
	uint64_t j = 0;
	uint64_t common = 0;
	uint64_t denom = 0;
	const HashList & hashesSortedRef = this->reference.hashesSorted;
	const HashList & hashesSortedQry = msh->reference.hashesSorted;
	//cout << "the size of hashesSortedRef is: " << this->reference.hashesSorted.size() << endl;
	//cout << "the size of hashesSortedQry is: " << msh->reference.hashesSorted.size() << endl;

	while ( denom < sketchSize && i < hashesSortedRef.size() && j < hashesSortedQry.size() )
	{
		if ( hashLessThan(hashesSortedRef.at(i), hashesSortedQry.at(j), hashesSortedRef.get64()) )
		{
			i++;
		}
		else if ( hashLessThan(hashesSortedQry.at(j), hashesSortedRef.at(i), hashesSortedRef.get64()) )
		{
			j++;
		}
		else
		{
			i++;
			j++;
			common++;
		}

		denom++;
	}

	if ( denom < sketchSize )
	{
		// complete the union operation if possible

		if ( i < hashesSortedRef.size() )
		{
			denom += hashesSortedRef.size() - i;
		}

		if ( j < hashesSortedQry.size() )
		{
			denom += hashesSortedQry.size() - j;
		}

		if ( denom > sketchSize )
		{
			denom = sketchSize;
		}
	}

	//	cout << "the common is: " << common << endl;
	//	cout << "the denom is: " << denom << endl;

	double jaccard = double(common) / denom;
	return jaccard;


}

double MinHash::dist(MinHash * msh)
{
	double distance;
	double maxDistance = 1;
	double maxPValue = 1;

	uint64_t sketchSize = this->parameters.minHashesPerWindow;
	int kmerSize = this->parameters.kmerSize;
	//int kmerSpace = this->kmerSpace;

	uint64_t i = 0;
	uint64_t j = 0;
	uint64_t common = 0;
	uint64_t denom = 0;
	const HashList & hashesSortedRef = this->reference.hashesSorted;
	const HashList & hashesSortedQry = msh->reference.hashesSorted;
	//	cout << "the size of hashesSortedRef is: " << hashesSortedRef.size() << endl;
	//	cout << "the size of hashesSortedQry is: " << hashesSortedQry.size() << endl;
	//	cout << "the size of hashesSortedRef is: " << this->reference.hashesSorted.size() << endl;
	//	cout << "the size of hashesSortedQry is: " << msh->reference.hashesSorted.size() << endl;

	while ( denom < sketchSize && i < hashesSortedRef.size() && j < hashesSortedQry.size() )
	{
		if ( hashLessThan(hashesSortedRef.at(i), hashesSortedQry.at(j), hashesSortedRef.get64()) )
		{
			i++;
		}
		else if ( hashLessThan(hashesSortedQry.at(j), hashesSortedRef.at(i), hashesSortedRef.get64()) )
		{
			j++;
		}
		else
		{
			i++;
			j++;
			common++;
		}

		denom++;
	}

	if ( denom < sketchSize )
	{
		// complete the union operation if possible

		if ( i < hashesSortedRef.size() )
		{
			denom += hashesSortedRef.size() - i;
		}

		if ( j < hashesSortedQry.size() )
		{
			denom += hashesSortedQry.size() - j;
		}

		if ( denom > sketchSize )
		{
			denom = sketchSize;
		}
	}



	double jaccard_ = double(common) / denom;
	distance = -log(2 * jaccard_ / (1. + jaccard_)) / kmerSize;

	if ( distance > 1 )
	{
		distance = 1;
	}

	if ( maxDistance >= 0 && distance > maxDistance )
	{
		return 1.;
	}
	double pValue_ = pValue(common, this->length, msh->length, kmerSpace, denom);
	if ( maxPValue >= 0 && pValue_ > maxPValue )
	{
		cerr << "the pValue is larger than maxPValue " << endl;
		return 1.;
	}
	return distance;


}

double MinHash::pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize)
{
	if ( x == 0 )
	{
		return 1.;
	}

	double pX = 1. / (1. + kmerSpace / lengthRef);
	double pY = 1. / (1. + kmerSpace / lengthQuery);
	//  cerr << endl;
	//	cerr << "kmerspace: " << kmerSpace << endl;
	//	cerr << "px: " << pX << endl;
	//	cerr << "py: " << pY << endl;
	double r = pX * pY / (pX + pY - pX * pY);

	//double M = (double)kmerSpace * (pX + pY) / (1. + r);

	//return gsl_cdf_hypergeometric_Q(x - 1, r * M, M - r * M, sketchSize);
	// 	return 0.1;   
	//	cerr << "r = " << r << endl;
	return gsl_cdf_binomial_Q(x - 1, r, sketchSize);
	//#ifdef USE_BOOST
	//    return cdf(complement(binomial(sketchSize, r), x - 1));
	//#else
	//    return gsl_cdf_binomial_Q(x - 1, r, sketchSize);
	//#endif
}
//start the WMinHash@xxm
WMinHash::WMinHash(Parameters parametersNew):parameters(parametersNew)
{
	binsArr = (double *) malloc (parameters.numBins * sizeof(double));
	for(int i = 0; i < parameters.numBins; i++){
		binsArr[i] = 0.0;
	}

	int g = ceil(2 / EPSILON);
	int d = ceil(log(1 - DELTA) / log(0.5));
	countMinSketch = (double *) malloc (d * g *sizeof(double));
	for(int i = 0; i < d * g; i++){
		countMinSketch[i] = 0.0;
	}

	r = (double *) malloc (parameters.histoSketch_sketchSize * parameters.histoSketch_dimension * sizeof(double));
	c = (double *) malloc (parameters.histoSketch_sketchSize * parameters.histoSketch_dimension * sizeof(double));
	b = (double *) malloc (parameters.histoSketch_sketchSize * parameters.histoSketch_dimension * sizeof(double));
	getCWS(r, c, b, parameters.histoSketch_sketchSize, parameters.histoSketch_dimension);
	
	histoSketch_sketch = (uint32_t *) malloc (parameters.histoSketch_sketchSize * sizeof(uint32_t));
	histoSketch_sketchWeight = (double *) malloc (parameters.histoSketch_sketchSize * sizeof(double));

	needToCompute = true;

}

void WMinHash::update(char * seq)
{
	int k = parameters.kmerSize;
	int w = parameters.minimizerWindowSize;
	int numBins = parameters.numBins;
	int histoSketch_sketchSize = parameters.histoSketch_sketchSize;
	int histoSketch_dimension = parameters.histoSketch_dimension;

	findMinimizers(k, w, seq, sketches);
	kmerSpectrumAddHash(sketches, binsArr, parameters.numBins);
//	for(int i = 0; i < parameters.numBins; i++){
//		if(binsArr[i] != 0){
//			printf("%d\t%lf\n", i, binsArr[i]);
//		}
//	}
	sketches.clear();
	needToCompute = true;

}

void WMinHash::computeHistoSketch()
{
	kmerSpectrumDump(binsArr, parameters.numBins, kmerSpectrums);
	for(int i = 0; i < kmerSpectrums.size(); i++){
		histoSketchAddElement((uint64_t)kmerSpectrums[i].BinID, kmerSpectrums[i].Frequency, countMinSketch, parameters.histoSketch_sketchSize, false, r, c, b, parameters.histoSketch_sketchSize, parameters.histoSketch_dimension, histoSketch_sketch, histoSketch_sketchWeight);
	}

}

void WMinHash::getWMinHash(){
	if(needToCompute){
		computeHistoSketch();
		needToCompute = false;
	}
	cout << "the sketch is: " << endl;
	for(int i = 0; i < parameters.histoSketch_sketchSize; i++){
		printf("%d\n", histoSketch_sketch[i]);
	}

	cout << "the sketchweight is: " << endl;
	for(int i = 0; i < parameters.histoSketch_sketchSize; i++){
		printf("%lf\n", histoSketch_sketchWeight[i]);
	}

}

double WMinHash::wJaccard(WMinHash * wmh)
{
	cout << "the needToCompute is: " << needToCompute << endl;
	if(needToCompute){
		computeHistoSketch();
		needToCompute = false;
	}
	if(wmh->needToCompute){
		wmh->computeHistoSketch();
		wmh->needToCompute = false;
	}
	double intersectElement = 0.0;
	double unionElement = 0.0;
	double jaccard = 0.0;
	for(int i = 0 ; i < parameters.histoSketch_sketchSize; i++){
		double curWeightA = abs(histoSketch_sketchWeight[i]);
		double curWeightB = abs(wmh->histoSketch_sketchWeight[i]);

		//get the intersection and union values
		if(histoSketch_sketch[i] == wmh->histoSketch_sketch[i]){
			if(curWeightA < curWeightB){
				intersectElement += curWeightA;
				unionElement += curWeightB;
			}
			else{
				intersectElement += curWeightB;
				unionElement += curWeightA;
			}
		}
		else{
			if(curWeightA > curWeightB){
				unionElement += curWeightA;
			}
			else{
				unionElement += curWeightB;
			}
		}
	
	}
	jaccard = intersectElement / unionElement;

	return jaccard;

}

	
double WMinHash::distance(WMinHash * wmh){
	return 1 - wJaccard(wmh);
}

WMinHash::~WMinHash()
{
	free(binsArr);
	free(countMinSketch);
	free(r);
	free(c);
	free(b);
	free(histoSketch_sketch);
	free(histoSketch_sketchWeight);

}

	


	








OMinHash::OMinHash(Parameters parametersNew, char * seqNew):
	parameters(parametersNew),
	seq(seqNew)
{
	m_k = parameters.kmerSize;
	m_l = parameters.l;
	m_m = parameters.m;
	rc = parameters.rc;

	if(rc){
		//TODO:get reverse complement to rcseq
	}
	sketch();
}

void OMinHash::sketch(){
	sk.k = m_k;		
	sk.l = m_l;		
	sk.m = m_m;		

	sk.data.resize(std::max(sk.l, 1) * sk.m * sk.k);
	compute_sketch(sk.data.data(), this->seq);

	if(rc){
		sk.rcdata.resize(std::max(sk.l, 1) * sk.m * sk.k);
		compute_sketch(sk.rcdata.data(), this->rcseq);
	}
}

inline void OMinHash::compute_sketch(char * ptr, const char * seq){
	std::string seqStr(seq);
	omh_pos(seqStr, m_k, m_l, m_m,
			[&ptr, &seq, this](unsigned i, unsigned j, size_t pos) { memcpy(ptr, seq + pos, m_k); ptr += m_k; });
}

template<typename BT>
static void omh_pos(const std::string& seq, unsigned k, unsigned l, unsigned m, BT block) {
	if(seq.size() < k) return;
	const bool weight = l > 0;
	if(l == 0) l = 1;

	std::vector<mer_info> mers;
	std::unordered_map<std::string, unsigned> occurrences;
	size_t pos[l];

	//  Create list of k-mers with occurrence numbers
	for(size_t i = 0; i < seq.size() - k + 1; ++i) {
		auto occ = occurrences[seq.substr(i, k)]++;
		mers.emplace_back(i, occ, (uint64_t)0);
	}

	xxhash hash;
	std::mt19937 gen32(32); //TODO: make 32 a parameter
	for(unsigned i = 0; i < m; ++i) {
		const auto seed = gen32();//prg();
		for(auto& meri : mers) {
			hash.reset(seed);
			hash.update(&seq.data()[meri.pos], k);
			if(weight) hash.update(&meri.occ, sizeof(meri.occ));
			meri.hash = hash.digest();
		}

		std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; });
		std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });
		for(unsigned j = 0; j < l; ++j)
			block(i, j, mers[j].pos);
	}
}

double OMinHash::compare_sketches(const OSketch& sk1, const OSketch& sk2, ssize_t m, bool circular) {
	if(sk1.k != sk2.k || sk1.l != sk2.l) return -1; // Different k or l
	if(m < 0) m = std::min(sk1.m, sk2.m);
	if(m > sk1.m || m > sk2.m) return -1;  // Too short

	const unsigned block = std::max(sk1.l, 1) * sk1.k;
	if(sk1.data.size() < m * block || sk2.data.size() < m * block) return -1; // Truncated
	const double fwd_score = compare_sketch_pair(sk1.data.data(), sk2.data.data(), m, sk1.k, sk1.l, circular);

	double bwd_score = 0.0;
	if(!sk1.rcdata.empty()) {
		bwd_score = compare_sketch_pair(sk1.rcdata.data(), sk2.data.data(), m, sk1.k, sk1.l, circular);
	} else if(!sk2.rcdata.empty()) {
		bwd_score = compare_sketch_pair(sk1.data.data(), sk2.rcdata.data(), m, sk1.k, sk1.l, circular);
	}

	return std::max(fwd_score, bwd_score);
}

double OMinHash::compare_sketch_pair(const char* p1, const char* p2, unsigned m, unsigned k, unsigned l, bool circular) {
	const unsigned block = std::max(l, (unsigned)1) * k;
	unsigned count = 0;
	if(!circular || l < 2) {
		for(unsigned i = 0; i < m; ++i, p1 += block, p2 += block)
			count += memcmp(p1, p2, block) == 0;
	} else {
		for(unsigned i = 0; i < m; ++i, p1 += block, p2 += block) {
			for(unsigned j = 0; j < l; ++j) {
				if(memcmp(p1, p2 + j * k, block - j * k) == 0 && memcmp(p1 + block - j * k, p2, j * k) == 0) {
					++count;
					break;
				}
			}
		}
	}
	return (double)count / m;
}

double OMinHash::similarity(OMinHash & omh2){
	return compare_sketches(this->sk, omh2.getSektch());	
}
