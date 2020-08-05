#include "MinHash.h"
#include "Sketch.h"
#include "MurmurHash3.h"
#include "hash.h"

#include <algorithm>
#include <iostream>
#include <immintrin.h>
#include <string.h>
#include <list>
#include <cstring>
#include <stdint.h>
#include <math.h>

#ifndef NOPYTHON
#include "pybind.h"
#endif

using namespace std;

hash_u HashList::at(int index) const
{
    hash_u hash;
    
    if ( use64 )
    {
        hash.hash64 = hashes64.at(index);
    }
    else
    {
        hash.hash32 = hashes32.at(index);
    }
    
    return hash;
}

void HashList::clear()
{
    if ( use64 )
    {
        hashes64.clear();
    }
    else
    {
        hashes32.clear();
    }
}

void HashList::resize(int size)
{
    if ( use64 )
    {
        hashes64.resize(size);
    }
    else
    {
        hashes32.resize(size);
    }
}

void HashList::set32(int index, uint32_t value)
{
    hashes32[index] = value;
}

void HashList::set64(int index, uint64_t value)
{
    hashes64[index] = value;
}

void HashList::sort()
{
    if ( use64 )
    {
        std::sort(hashes64.begin(), hashes64.end());
    }
    else
    {
        std::sort(hashes32.begin(), hashes32.end());
    }
}

void HashPriorityQueue::clear()
{
    if ( use64 )
    {
        while ( queue64.size() )
        {
            queue64.pop();
        }
    }
    else
    {
        while ( queue32.size() )
        {
            queue32.pop();
        }
    }
}

hash_u HashPriorityQueue::top() const
{
    hash_u hash;
    
    if ( use64 )
    {
        hash.hash64 = queue64.top();
    }
    else
    {
        hash.hash32 = queue32.top();
    }
    
    return hash;
}

uint32_t HashSet::count(hash_u hash) const
{
	if ( use64 )
	{
		if ( hashes64.count(hash.hash64) )
		{
			return hashes64.at(hash.hash64);
		}
		else
		{
			return 0;
		}
	}
	else
	{
		if ( hashes32.count(hash.hash32) )
		{
			return hashes32.at(hash.hash32);
		}
		else
		{
			return 0;
		}
	}
}

void HashSet::erase(hash_u hash)
{
    if ( use64 )
    {
        hashes64.erase(hash.hash64);
    }
    else
    {
        hashes32.erase(hash.hash32);
    }
}

void HashSet::insert(hash_u hash, uint32_t count)
{
    if ( use64 )
    {
    	hash64_t hash64 = hash.hash64;
    	
    	if ( hashes64.count(hash64) )
    	{
    		hashes64[hash64] = hashes64.at(hash64) + count;
    	}
    	else
    	{
	        hashes64[hash64] = count;
	    }
    }
    else
    {
    	hash32_t hash32 = hash.hash32;
    	
    	if ( hashes32.count(hash32) )
    	{
    		hashes32[hash32] = hashes32.at(hash32) + count;
    	}
    	else
    	{
	        hashes32[hash32] = count;
	    }
    }
}

void HashSet::toCounts(std::vector<uint32_t> & counts) const
{
    if ( use64 )
    {
        for ( robin_hood::unordered_map<hash64_t, uint32_t>::const_iterator i = hashes64.begin(); i != hashes64.end(); ++i )
        {
            counts.push_back(i->second);
        }
    }
    else
    {
        for ( robin_hood::unordered_map<hash32_t, uint32_t>::const_iterator i = hashes32.begin(); i != hashes32.end(); ++i )
        {
            counts.push_back(i->second);
        }
    }
}

void HashSet::toHashList(HashList & hashList) const
{
    if ( use64 )
    {
        for ( robin_hood::unordered_map<hash64_t, uint32_t>::const_iterator i = hashes64.begin(); i != hashes64.end(); ++i )
        {
            hashList.push_back64(i->first);
        }
    }
    else
    {
        for ( robin_hood::unordered_map<hash32_t, uint32_t>::const_iterator i = hashes32.begin(); i != hashes32.end(); ++i )
        {
            hashList.push_back32(i->first);
        }
    }
}

MinHashHeap::MinHashHeap(bool use64New, uint64_t cardinalityMaximumNew, uint64_t multiplicityMinimumNew, uint64_t memoryBoundBytes) :
	use64(use64New),
	hashes(use64New),
	hashesQueue(use64New),
	hashesPending(use64New),
	hashesQueuePending(use64New)
{
	cardinalityMaximum = cardinalityMaximumNew;
	multiplicityMinimum = multiplicityMinimumNew;
	
	multiplicitySum = 0;
	
}

MinHashHeap::~MinHashHeap()
{}

void MinHashHeap::computeStats()
{
	vector<uint32_t> counts;
	hashes.toCounts(counts);
	
	for ( int i = 0; i < counts.size(); i++ )
	{
		cout << counts.at(i) << endl;
	}
}

void MinHashHeap::clear()
{
	hashes.clear();
	hashesQueue.clear();
	
	hashesPending.clear();
	hashesQueuePending.clear();
	
	multiplicitySum = 0;
}

void MinHashHeap::tryInsert(hash_u hash)
{
	if
	(
		hashes.size() < cardinalityMaximum ||
		hashLessThan(hash, hashesQueue.top(), use64)
	)
	{
		if ( hashes.count(hash) == 0 )
		{

			if ( multiplicityMinimum == 1 || hashesPending.count(hash) == multiplicityMinimum - 1 )
			{
				hashes.insert(hash, multiplicityMinimum);
				hashesQueue.push(hash);
				multiplicitySum += multiplicityMinimum;
				
				if ( multiplicityMinimum > 1 )
				{
					// just remove from set for now; will be removed from
					// priority queue when it's on top
					//
					hashesPending.erase(hash);
				}
			}
			else
			{
				if ( hashesPending.count(hash) == 0 )
				{
					hashesQueuePending.push(hash);
				}
			
				hashesPending.insert(hash, 1);
			}
		}
		else
		{
			hashes.insert(hash, 1);
			multiplicitySum++;
		}
		
		if ( hashes.size() > cardinalityMaximum )
		{
			multiplicitySum -= hashes.count(hashesQueue.top());
			hashes.erase(hashesQueue.top());
			
			// loop since there could be zombie hashes (gone from hashesPending)
			//
			while ( hashesQueuePending.size() > 0 && hashLessThan(hashesQueue.top(), hashesQueuePending.top(), use64) )
			{
				if ( hashesPending.count(hashesQueuePending.top()) )
				{
					hashesPending.erase(hashesQueuePending.top());
				}
				
				hashesQueuePending.pop();
			}
			
			hashesQueue.pop();
		}
	}
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


namespace Sketch
{

#if defined __AVX512F__ && defined __AVX512CD__
__m512i inline min512(__m512i v1, __m512i v2){
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
#else 
#ifdef __AVX2__
	// implement by avx2
void inline transpose4_epi64(__m256i *row1, __m256i *row2, __m256i *row3, __m256i *row4)
{
	__m256i vt1, vt2, vt3, vt4;

	vt1 = _mm256_unpacklo_epi64(*row1, *row2);
	vt2 = _mm256_unpackhi_epi64(*row1, *row2);
	vt3 = _mm256_unpacklo_epi64(*row3, *row4);
	vt4 = _mm256_unpackhi_epi64(*row3, *row4);

	*row1 =_mm256_permute2x128_si256(vt1, vt3, 0x20);
	*row2 =_mm256_permute2x128_si256(vt2, vt4, 0x20);
	*row3 =_mm256_permute2x128_si256(vt1, vt3, 0x31);
	*row4 =_mm256_permute2x128_si256(vt2, vt4, 0x31);
}

	#else
	//non vectorization 
	#endif
#endif


//MinHash::MinHash()
//{
//	minHashHeap = new MinHashHeap(use64, sketchSize);// parameters.reads ?  parameters.minCov : 1);
//
//	this->kmerSpace = pow(alphabetSize, kmerSize);
//	//cerr << "kmerSpace init from pow is " << this->kmerSpace << endl;
//	this->totalLength = 0;
//	this->needToList = true;
//
//}

void MinHash::update(char * seq)
{
	const uint64_t length = strlen(seq);
	totalLength += length;

	// to uppercase 
	// note: kmers sequences with different cases will lead to different hashes
	if( ! preserveCase )
		for ( uint64_t i = 0; i < length; i++ )
		{
			if ( seq[i] > 96 && seq[i] < 123 )
			{
				seq[i] -= 32;
			}
		}

   char * seqRev;
    
    if ( ! noncanonical )
    {
    	seqRev = new char[length];
        reverseComplement(seq, seqRev, length);
    }

#if defined __AVX512F__ && defined __AVX512BW__

	int pend_k = ((kmerSize - 1) / 16 + 1) * 16;
	int n_kmers = length - kmerSize + 1;
	int n_kmers_body = (n_kmers / 16) * 16;
	
	const uint8_t* input8 = (const uint8_t *)seq;
	const uint8_t* input8_rev = (const uint8_t *)seqRev;
	uint64_t res[16 * 2];
	uint64_t res2[2];
	uint8_t kmer_buf[kmerSize];
	
	__m512i v0, v1;
	__m512i vi[8];
	__m512i vj[8];
	__m512i vi_forword[8];
	__m512i vi_reverse[8];
	__m512i vj_forword[8];
	__m512i vj_reverse[8];
	__m512i vzero = _mm512_setzero_si512();

	__mmask64 mask_load = 0xffffffffffffffff;
	mask_load >>= (64 - kmerSize);

	for(int i = 0; i < n_kmers_body-1; i+=16)
	{
		vi_forword[0] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 0);
		vi_forword[1] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 1);
		vi_forword[2] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 2);
		vi_forword[3] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 3);
		vi_forword[4] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 4);
		vi_forword[5] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 5);
		vi_forword[6] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 6);
		vi_forword[7] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 7);


		vj_forword[0] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 8);
		vj_forword[1] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 9);
		vj_forword[2] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 10);
		vj_forword[3] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 11);
		vj_forword[4] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 12);
		vj_forword[5] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 13);
		vj_forword[6] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 14);
		vj_forword[7] = _mm512_mask_loadu_epi8(vzero, mask_load, input8 + i + 15);

		if( !noncanonical){

			vi_reverse[0] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 0 - kmerSize);
			vi_reverse[1] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 1 - kmerSize);
			vi_reverse[2] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 2 - kmerSize);
			vi_reverse[3] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 3 - kmerSize);
			vi_reverse[4] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 4 - kmerSize);
			vi_reverse[5] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 5 - kmerSize);
			vi_reverse[6] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 6 - kmerSize);
			vi_reverse[7] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 7 - kmerSize);

			vj_reverse[0] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 8 - kmerSize);
			vj_reverse[1] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 9 - kmerSize);
			vj_reverse[2] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 10 - kmerSize);
			vj_reverse[3] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 11 - kmerSize);
			vj_reverse[4] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 12 - kmerSize);
			vj_reverse[5] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 13 - kmerSize);
			vj_reverse[6] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 14 - kmerSize);
			vj_reverse[7] = _mm512_mask_loadu_epi8(vzero, mask_load, input8_rev + length - i - 15 - kmerSize);


			vi[0] = min512(vi_forword[0], vi_reverse[0]);
			vi[1] = min512(vi_forword[1], vi_reverse[1]);
			vi[2] = min512(vi_forword[2], vi_reverse[2]);
			vi[3] = min512(vi_forword[3], vi_reverse[3]);
			vi[4] = min512(vi_forword[4], vi_reverse[4]);
			vi[5] = min512(vi_forword[5], vi_reverse[5]);
			vi[6] = min512(vi_forword[6], vi_reverse[6]);
			vi[7] = min512(vi_forword[7], vi_reverse[7]);

			vj[0] = min512(vj_forword[0], vj_reverse[0]);
			vj[1] = min512(vj_forword[1], vj_reverse[1]);
			vj[2] = min512(vj_forword[2], vj_reverse[2]);
			vj[3] = min512(vj_forword[3], vj_reverse[3]);
			vj[4] = min512(vj_forword[4], vj_reverse[4]);
			vj[5] = min512(vj_forword[5], vj_reverse[5]);
			vj[6] = min512(vj_forword[6], vj_reverse[6]);
			vj[7] = min512(vj_forword[7], vj_reverse[7]);

		}else{

			vi[0] = vi_forword[0];
			vi[1] = vi_forword[1];
			vi[2] = vi_forword[2];
			vi[3] = vi_forword[3];
			vi[4] = vi_forword[4];
			vi[5] = vi_forword[5];
			vi[6] = vi_forword[6];
			vi[7] = vi_forword[7];

			vj[0] = vj_forword[0];
			vj[1] = vj_forword[1];
			vj[2] = vj_forword[2];
			vj[3] = vj_forword[3];
			vj[4] = vj_forword[4];
			vj[5] = vj_forword[5];
			vj[6] = vj_forword[6];
			vj[7] = vj_forword[7];
		}

	//	vi[0] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 0, input8_rev + length - i - 0 - kmerSize, kmerSize) <= 0 ? input8 + i + 0 : input8_rev + length - i - 0 - kmerSize);
	//	vi[1] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 1, input8_rev + length - i - 1 - kmerSize, kmerSize) <= 0 ? input8 + i + 1 : input8_rev + length - i - 1 - kmerSize);
	//	vi[2] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 2, input8_rev + length - i - 2 - kmerSize, kmerSize) <= 0 ? input8 + i + 2 : input8_rev + length - i - 2 - kmerSize);
	//	vi[3] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 3, input8_rev + length - i - 3 - kmerSize, kmerSize) <= 0 ? input8 + i + 3 : input8_rev + length - i - 3 - kmerSize);
	//	vi[4] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 4, input8_rev + length - i - 4 - kmerSize, kmerSize) <= 0 ? input8 + i + 4 : input8_rev + length - i - 4 - kmerSize);
	//	vi[5] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 5, input8_rev + length - i - 5 - kmerSize, kmerSize) <= 0 ? input8 + i + 5 : input8_rev + length - i - 5 - kmerSize);
	//	vi[6] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 6, input8_rev + length - i - 6 - kmerSize, kmerSize) <= 0 ? input8 + i + 6 : input8_rev + length - i - 6 - kmerSize);
	//	vi[7] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 7, input8_rev + length - i - 7 - kmerSize, kmerSize) <= 0 ? input8 + i + 7 : input8_rev + length - i - 7 - kmerSize);

	//	vj[0] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 8, input8_rev + length - i - 8 - kmerSize, kmerSize) <= 0 ? input8 + i + 8 : input8_rev + length - i - 8 - kmerSize);
	//	vj[1] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 9, input8_rev + length - i - 9 - kmerSize, kmerSize) <= 0 ? input8 + i + 9 : input8_rev + length - i - 9 - kmerSize);
	//	vj[2] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 10, input8_rev + length - i - 10 - kmerSize, kmerSize) <= 0 ? input8 + i + 10 : input8_rev + length - i - 10 - kmerSize);
	//	vj[3] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 11, input8_rev + length - i - 11 - kmerSize, kmerSize) <= 0 ? input8 + i + 11 : input8_rev + length - i - 11 - kmerSize);
	//	vj[4] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 12, input8_rev + length - i - 12 - kmerSize, kmerSize) <= 0 ? input8 + i + 12 : input8_rev + length - i - 12 - kmerSize);
	//	vj[5] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 13, input8_rev + length - i - 13 - kmerSize, kmerSize) <= 0 ? input8 + i + 13 : input8_rev + length - i - 13 - kmerSize);
	//	vj[6] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 14, input8_rev + length - i - 14 - kmerSize, kmerSize) <= 0 ? input8 + i + 14 : input8_rev + length - i - 14 - kmerSize);
	//	vj[7] = _mm512_mask_loadu_epi8(vzero, mask_load, memcmp(input8 + i + 15, input8_rev + length - i - 15 - kmerSize, kmerSize) <= 0 ? input8 + i + 15 : input8_rev + length - i - 15 - kmerSize);



		transpose8_epi64(&vi[0], &vi[1], &vi[2], &vi[3], &vi[4], &vi[5], &vi[6], &vi[7]); 
		transpose8_epi64(&vj[0], &vj[1], &vj[2], &vj[3], &vj[4], &vj[5], &vj[6], &vj[7]); 

		MurmurHash3_x64_128_avx512_8x16(vi, vj, pend_k, kmerSize, seed, res);// the seed in Mash is 42; verified by xxm;

		hash_u hash;
		for(int j = 0; j < 16; j++){
			if(use64)
				hash.hash64 = res[j * 2];
			else
				hash.hash32 = (uint32_t)res[j * 2];
			minHashHeap->tryInsert(hash);
		}
	}

	//tail
	for(int i = n_kmers_body; i < n_kmers; i++){
		bool noRev = true;
		if(!noncanonical && (memcmp(input8 + i, input8_rev + length - i - kmerSize, kmerSize) >= 0)){
			noRev = false;
		}
		if(noRev)
		{
			for(int j = 0; j < kmerSize; j++)
			{
				//kmer_buf[j] = input8[i + j];
				kmer_buf[j] = input8[i + j];
			}
		}
		else
		{
			for(int j = 0; j < kmerSize; j++)
			{
				//kmer_buf[j] = input8[i + j];
				kmer_buf[j] = input8_rev[length - i - kmerSize + j];
			}
		}

		MurmurHash3_x64_128(kmer_buf, kmerSize, seed, res2);// the getHash just need the lower 64bit of the total 128bit of res[i];
		hash_u hash;
		if(use64)
			hash.hash64 = res2[0];
		else
			hash.hash32 = (uint32_t)res2[0];

		minHashHeap->tryInsert(hash);

	}
//============================================================================================================================================================

#else
	#if defined __AVX2__

	int pend_k = ((kmerSize - 1) / 16 + 1) * 16;
	int n_kmers = length - kmerSize + 1;
	int n_kmers_body = (n_kmers / 4) * 4;
	const uint8_t * input8 = (const uint8_t *)seq;
	const uint8_t * input8_rev = (const uint8_t *)seqRev;
	uint64_t res[4 * 2];
	uint64_t res2[2];
	uint8_t kmer_buf[kmerSize];

	__m256i vi[4];
	//__m256i vzero = _mm256_set1_epi32(0x0);
	uint8_t maskArr[32];
	for(int i = 0; i < kmerSize; i++){
		maskArr[i] = 0xff;
	}
	for(int i = kmerSize; i < 32; i++){
		maskArr[i] = 0x0;
	}
	__m256i vmask = _mm256_loadu_si256((__m256i *)maskArr);
	for(int i = 0; i < n_kmers_body - 1; i+=4){
		vi[0] = _mm256_loadu_si256((__m256i *)(noncanonical || memcmp(input8 + i + 0, input8_rev + length - i - 0 - kmerSize, kmerSize) <= 0 ? input8 + i + 0 : input8_rev + length - i - 0 - kmerSize));
		vi[1] = _mm256_loadu_si256((__m256i *)(noncanonical || memcmp(input8 + i + 1, input8_rev + length - i - 1 - kmerSize, kmerSize) <= 0 ? input8 + i + 1 : input8_rev + length - i - 1 - kmerSize));
		vi[2] = _mm256_loadu_si256((__m256i *)(noncanonical || memcmp(input8 + i + 2, input8_rev + length - i - 2 - kmerSize, kmerSize) <= 0 ? input8 + i + 2 : input8_rev + length - i - 2 - kmerSize));
		vi[3] = _mm256_loadu_si256((__m256i *)(noncanonical || memcmp(input8 + i + 3, input8_rev + length - i - 3 - kmerSize, kmerSize) <= 0 ? input8 + i + 3 : input8_rev + length - i - 3 - kmerSize));
		vi[0] = _mm256_and_si256(vi[0], vmask);
		vi[1] = _mm256_and_si256(vi[1], vmask);
		vi[2] = _mm256_and_si256(vi[2], vmask);
		vi[3] = _mm256_and_si256(vi[3], vmask);

		transpose4_epi64(&vi[0], &vi[1], &vi[2], &vi[3]);

		MurmurHash3_x64_128_avx2_8x4(vi, pend_k, kmerSize, seed, res);

		hash_u hash;
		for(int j = 0; j < 4; j++)
		{
			if(use64)	
				hash.hash64 = res[j * 2];
			else
				hash.hash32 = (uint32_t)res[j * 2];

			minHashHeap->tryInsert(hash);
		}
	}

	for(int i = n_kmers_body; i < n_kmers; i++)
	{
		//bool noRev = (memcmp(input8 + i, input8_rev + length -i - kmerSize, kmerSize) <= 0) || noncanonical;
		bool noRev = true;
		if(!noncanonical && (memcmp(input8 + i, input8_rev + length - i - kmerSize, kmerSize) >= 0)){
			noRev = false;
		}
		if(noRev){
			for(int j = 0; j < kmerSize; j++){
				kmer_buf[j] = input8[i + j];
			}
		}
		else{
			for(int j = 0; j < kmerSize; j++){
				kmer_buf[j] = input8_rev[length - i - kmerSize + j];
			}
		}

		MurmurHash3_x64_128(kmer_buf, kmerSize, seed, res2);
		hash_u hash;
		if(use64)
			hash.hash64 = res2[0];
		else
			hash.hash32 = (uint32_t)res2[0];
			
		minHashHeap->tryInsert(hash);

	}

	#else

	//implement by no optmization
    for ( uint64_t i = 0; i < length - kmerSize + 1; i++ )
    {
            
        const char *kmer_fwd = seq + i;
        const char *kmer_rev = seqRev + length - i - kmerSize;
        const char * kmer = (noncanonical || memcmp(kmer_fwd, kmer_rev, kmerSize) <= 0) ? kmer_fwd : kmer_rev;
        bool filter = false;
        
        hash_u hash = getHash(kmer, kmerSize, seed, use64);
        
		minHashHeap->tryInsert(hash);
    }
	#endif
#endif
    
    
    if ( ! noncanonical )
    {
        delete [] seqRev;
    }


	needToList = true;
}

void MinHash::heapToList()
{
	HashList & hashlist = reference.hashesSorted;
	hashlist.clear();
	hashlist.setUse64(use64);
	minHashHeap -> toHashList(hashlist);
	minHashHeap -> toCounts(reference.counts);
	hashlist.sort();

}

void MinHash::printMinHashes()
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
		if(use64)
			cerr << "hash64 " <<  i << " " << reference.hashesSorted.at(i).hash64 << endl;
		else
			cerr << "hash32 " <<  i << " " << reference.hashesSorted.at(i).hash32 << endl;
	}
	return;
}


void MinHash::merge(MinHash& msh)
{
	msh.heapToList();
	HashList & mshList = msh.reference.hashesSorted;	
	for(int i = 0; i < mshList.size(); i++)
	{
		//cerr << "insert to heap" << mshList.at(i).hash64 << endl;
		minHashHeap -> tryInsert(mshList.at(i));
	}
	needToList = true;
		
	return;	
}

double MinHash::jaccard(MinHash * msh)
{
	if(needToList){
		heapToList();
		needToList = false;
	}
	if(msh->needToList){
		//cout << "msh2 need to list addbyxxm " << endl;
		msh->heapToList();
		msh->needToList = false;
	}

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

double MinHash::mdistance(MinHash * msh)
{
	double distance;
	double maxDistance = 1;
	double maxPValue = 1;

	double jaccard_ = this->jaccard(msh);
	distance = -log(2 * jaccard_ / (1. + jaccard_)) / kmerSize;

	if ( distance > 1 )
	{
		distance = 1;
	}

	if ( maxDistance >= 0 && distance > maxDistance )
	{
		return 1.;
	}
	//double pValue_ = pValue(common, this->length, msh->length, kmerSpace, denom);
	//if ( maxPValue >= 0 && pValue_ > maxPValue )
	//{
	//	cerr << "the pValue is larger than maxPValue " << endl;
	//	return 1.;
	//}
	return distance;


}
/*
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
*/

//FIXME: it only works for DNA sequences
void reverseComplement(const char * src, char * dest, int length)
{
	char table[4] = {'T','G','A','C'};
    for ( int i = 0; i < length; i++ )
    {
        char base = src[i];
		base >>= 1;
		base &= 0x03;
        dest[length - i - 1] = table[base];
    }
}

}
