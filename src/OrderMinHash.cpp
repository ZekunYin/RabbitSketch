
#include "OrderMinHash.h"
#include "Sketch.h"
#include "xxhash.hpp"
#include "hash_int.h"

#include <algorithm>
#include <queue>

#include "robin_hood.h"

#include <sys/time.h>

#include <immintrin.h>

double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}
#ifdef __AVX512F__
static void inspect(__m512i va)
{
	uint64_t a[8];
	_mm512_storeu_si512(a, va);
	cout << "[ ";
	for(int i = 0; i<8;i++)
		cout << a[i] << " ";
	cout << "]" << endl;
}
#endif
namespace Sketch
{

OrderMinHash::OrderMinHash(char * seqNew):
	seq(seqNew)
{

	if(rc){
		//TODO:get reverse complement to rcseq
		int rc_len = strlen(seq);
		rcseq = new char[rc_len];
		char table[4] = {'T','G','A','C'};
		for ( uint64_t i = 0; i < rc_len; i++ )
		{
			char base = seq[i];

			base >>= 1;
			base &= 0x03;
			rcseq[rc_len - i - 1] = table[base];

		}
	
	}
	sketch();
}

void OrderMinHash::buildSketch(char * seqNew = NULL)
{
	// rebuild sketch using old data
	if(seqNew == NULL)
	{
		if(seq == NULL)
		{
			cerr << "WARNING: no data found" << endl;
			return;
		}

		if(rc){
			int rc_len = strlen(seq);

			if(rcseq != NULL){
				delete [] rcseq;
				rcseq = NULL;
			}

			rcseq = new char[rc_len];
			reverseComplement(seq, rcseq, rc_len);
		}
		sketch();
	} else {
		seq = seqNew;

		if(rc){
			int rc_len = strlen(seq);

			if(rcseq != NULL){
				delete [] rcseq;
				rcseq = NULL;
			}

			rcseq = new char[rc_len];
			reverseComplement(seq, rcseq, rc_len);
		
		}
		sketch();

	}
}
void OrderMinHash::sketch()
{
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

inline void OrderMinHash::compute_sketch(char * ptr, const char * seq){
	std::string seqStr(seq);
	double t1 = get_sec();
	omh_pos(seqStr, m_k, m_l, m_m, mtSeed, ptr);
	double t2 = get_sec();
	cout << "total sketch time: " << t2 - t1 << endl;
//			[&ptr, &seq, this](unsigned i, unsigned j, size_t pos) { memcpy(ptr, seq + pos, m_k); ptr += m_k; });
}

//template<typename BT>
static void omh_pos(const std::string& seq, unsigned k, unsigned l, unsigned m, uint64_t mtSeed, char * ptr) {
	if(seq.size() < k) return;
	const uint64_t weight = l > 0 ? 1 : 0;
	if(l == 0) l = 1;

	std::vector<mer_info> mers;
	robin_hood::unordered_map<std::string, unsigned> occurrences;
	size_t *pos = new size_t[l];

	//  Create list of k-mers with occurrence numbers
	double t1 = get_sec();
	for(size_t i = 0; i < seq.size() - k + 1; ++i) {
		auto occ = occurrences[seq.substr(i, k)]++;
		//mers.emplace_back(i, occ, (uint64_t)0, hash_to_uint(&seq.data()[i], k));
		mers.emplace_back(i, occ, (uint64_t)0, 0);
	}
	uint64_t * intHash = new uint64_t[seq.size() -k + 1];
	uint64_t * occ     = new uint64_t[seq.size() -k + 1];
	for(int i = 0; i < seq.size() -k + 1; i++)
	{
		intHash[i] = hash_to_uint(&seq.data()[i], k);
		occ[i]     = mers[i].occ;
	}
	double t2 = get_sec();
	std::cout << "occ time: " << t2 - t1 << std::endl;
	std::cout << "omh mers size:" << mers.size() << std::endl;
	std::cout << "seq size: " << seq.size() << std::endl;
	std::cout << "seq -k +1: " << seq.size() - k + 1 << std::endl;
	xxhash hash;
	std::mt19937_64 gen64(mtSeed); //TODO: make 32 a parameter
	t1 = get_sec();
	auto cmp = [](Sketch::mer_info & a, Sketch::mer_info & b){return a.hash < b.hash;};
	std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> pqueue(cmp);
	std::vector<mer_info> lmers;
	lmers.reserve(l);
	int lanes = 8*4;// vector lanes unroll(4)
	uint64_t buffer[lanes];
	mer_info * minmers = new mer_info[l];

#ifdef __AVX512F__
	__m512i va[4];
	__m512i vb[4];
	__m512i vtmp[4];
	__m512i vocc[4];
	__m512i vweight = _mm512_set1_epi64(weight);
	__m512i vzero   = _mm512_set1_epi64(0);
	__m512i vconst0 = _mm512_set1_epi64(0xff51afd7ed558ccd);
	__m512i vconst1 = _mm512_set1_epi64(0xc4ceb9fe1a85ec53);
	__mmask8 weight_msk = l > 0 ? 0xFF : 0x00;
#endif

	for(unsigned i = 0; i < m; ++i) {
		const auto seed = gen64();//prg();
		for(int id = 0; id < l; id++)
			minmers[id].hash = UINT64_MAX;

/*
#ifdef __AVX512F__
		__m512i vseed = _mm512_set1_epi64(seed);
#endif

		size_t pend_size = (mers.size() / lanes) * lanes;
		//body
		for( int id = 0; id < pend_size; id+=lanes)
		{
			//uint64_t kmer_int = mers[id].int_hash;
#ifdef __NOAVX512F__ 
		va[0] = _mm512_loadu_si512((void *)&intHash[id + 0 * 8]);
		va[1] = _mm512_loadu_si512((void *)&intHash[id + 1 * 8]);
		va[2] = _mm512_loadu_si512((void *)&intHash[id + 2 * 8]);
		va[3] = _mm512_loadu_si512((void *)&intHash[id + 3 * 8]);

		//vocc[0] = _mm512_mask_loadu_epi64(vzero, weight_msk, (void *)&occ[id + 0 * 8]);
		//vocc[1] = _mm512_mask_loadu_epi64(vzero, weight_msk, (void *)&occ[id + 1 * 8]);
		//vocc[2] = _mm512_mask_loadu_epi64(vzero, weight_msk, (void *)&occ[id + 2 * 8]);
		//vocc[3] = _mm512_mask_loadu_epi64(vzero, weight_msk, (void *)&occ[id + 3 * 8]);


		vocc[0] = _mm512_loadu_si512((void *)&occ[id + 0 * 8]);
		vocc[1] = _mm512_loadu_si512((void *)&occ[id + 1 * 8]);
		vocc[2] = _mm512_loadu_si512((void *)&occ[id + 2 * 8]);
		vocc[3] = _mm512_loadu_si512((void *)&occ[id + 3 * 8]);

		vb[0] = _mm512_mask_add_epi64(va[0], weight_msk, va[0], vocc[0]);
		vb[1] = _mm512_mask_add_epi64(va[1], weight_msk, va[1], vocc[1]);
		vb[2] = _mm512_mask_add_epi64(va[2], weight_msk, va[2], vocc[2]);
		vb[3] = _mm512_mask_add_epi64(va[3], weight_msk, va[3], vocc[3]);

		//vtmp[0] = _mm512_mullo_epi64(vocc[0], vweight);
		//vtmp[1] = _mm512_mullo_epi64(vocc[1], vweight);
		//vtmp[2] = _mm512_mullo_epi64(vocc[2], vweight);
		//vtmp[3] = _mm512_mullo_epi64(vocc[3], vweight);

		//vb[0] = _mm512_add_epi64(va[0], vocc[0]);
		//vb[1] = _mm512_add_epi64(va[1], vocc[1]);
		//vb[2] = _mm512_add_epi64(va[2], vocc[2]);
		//vb[3] = _mm512_add_epi64(va[3], vocc[3]);
	

		va[0] = _mm512_xor_epi64(vb[0], vseed);
		va[1] = _mm512_xor_epi64(vb[1], vseed);
		va[2] = _mm512_xor_epi64(vb[2], vseed);
		va[3] = _mm512_xor_epi64(vb[3], vseed);

		vtmp[0] = _mm512_srli_epi64(va[0], 33);
		vtmp[1] = _mm512_srli_epi64(va[1], 33);
		vtmp[2] = _mm512_srli_epi64(va[2], 33);
		vtmp[3] = _mm512_srli_epi64(va[3], 33);

		vb[0] = _mm512_xor_epi64(va[0], vtmp[0]);
		vb[1] = _mm512_xor_epi64(va[1], vtmp[1]);
		vb[2] = _mm512_xor_epi64(va[2], vtmp[2]);
		vb[3] = _mm512_xor_epi64(va[3], vtmp[3]);

		va[0] = _mm512_mullo_epi64(vb[0], vconst0);
		va[1] = _mm512_mullo_epi64(vb[1], vconst0);
		va[2] = _mm512_mullo_epi64(vb[2], vconst0);
		va[3] = _mm512_mullo_epi64(vb[3], vconst0);

		vtmp[0] = _mm512_srli_epi64(va[0], 33);
		vtmp[1] = _mm512_srli_epi64(va[1], 33);
		vtmp[2] = _mm512_srli_epi64(va[2], 33);
		vtmp[3] = _mm512_srli_epi64(va[3], 33);

		vb[0] = _mm512_xor_epi64(va[0], vtmp[0]);
		vb[1] = _mm512_xor_epi64(va[1], vtmp[1]);
		vb[2] = _mm512_xor_epi64(va[2], vtmp[2]);
		vb[3] = _mm512_xor_epi64(va[3], vtmp[3]);

		va[0] = _mm512_mullo_epi64(vb[0], vconst1);
		va[1] = _mm512_mullo_epi64(vb[1], vconst1);
		va[2] = _mm512_mullo_epi64(vb[2], vconst1);
		va[3] = _mm512_mullo_epi64(vb[3], vconst1);

		vtmp[0] = _mm512_srli_epi64(va[0], 33);
		vtmp[1] = _mm512_srli_epi64(va[1], 33);
		vtmp[2] = _mm512_srli_epi64(va[2], 33);
		vtmp[3] = _mm512_srli_epi64(va[3], 33);

		vb[0] = _mm512_xor_epi64(va[0], vtmp[0]);
		vb[1] = _mm512_xor_epi64(va[1], vtmp[1]);
		vb[2] = _mm512_xor_epi64(va[2], vtmp[2]);
		vb[3] = _mm512_xor_epi64(va[3], vtmp[3]);

		_mm512_storeu_si512(buffer + 0 * 8, vb[0]);
		_mm512_storeu_si512(buffer + 1 * 8, vb[1]);
		_mm512_storeu_si512(buffer + 2 * 8, vb[2]);
		_mm512_storeu_si512(buffer + 3 * 8, vb[3]);
#else
			#pragma ivdep
			#pragma unroll(4)
			for(int vid = 0; vid < lanes; vid++)	
			{
				uint64_t kmer_int = intHash[id + vid];
		 	    kmer_int += occ[id + vid] * weight;
				buffer[vid] = kmer_int;
				buffer[vid] ^= seed; 
    			buffer[vid] ^= buffer[vid] >> 33;
    			buffer[vid] *= 0xff51afd7ed558ccd;
    			buffer[vid] ^= buffer[vid] >> 33;
    			buffer[vid] *= 0xc4ceb9fe1a85ec53;
    			buffer[vid] ^= buffer[vid] >> 33;

			}
#endif
			//meri.hash = hash.digest();
			for(int vid = 0; vid < lanes; vid++)	
			{
				//for(int index = 0; index < l; index ++)
				//{
				//	cout << minmers[index].hash << " ";
				//}
				//cout << endl;
				//cout << "insert hash: " << buffer[vid] << endl;

				mers[id + vid].hash = buffer[vid];
				//insert mers to minmers
				//if(mers[id+vid].hash < minmers[l - 1].hash)
				//{
				//	int break_id = 0;
				//	for(int insert_id = l-2; insert_id >= 0; insert_id--)
				//	{
				//		if(mers[id + vid].hash >= minmers[insert_id].hash) 
				//		{
				//			break_id = insert_id + 1;
				//			break;
				//		}	
				//	}
				//	//cout << "break id: " << break_id << endl;
				//	memcpy(minmers + break_id + 1, minmers + break_id, (l - break_id - 1) * sizeof(mer_info) );
				//	minmers[break_id] = mers[id + vid];
				//}
				//if(i == 0) cout << mers[id+vid].hash << endl; 
				//if(pqueue.empty())
				//	pqueue.push(mers[id + vid]);
				//else 
				if( pqueue.size() < l || mers[id + vid].hash < pqueue.top().hash)
				{
					pqueue.push(mers[id + vid]);
					if(pqueue.size() > l) pqueue.pop();
				}
			}
		}
		*/
		//exit(0);
		//tail
		for( int id = 0; id < mers.size(); id++)
		{

			uint64_t kmer_int = intHash[id];
		    kmer_int += occ[id] * weight;
			mers[id].hash = mc::murmur3_fmix(kmer_int, seed);
			//if(mers[id].hash < minmers[l - 1].hash)
			//{
			//	int break_id = 0;
			//	for(int insert_id = l-2; insert_id >= 0; insert_id--)
			//	{
			//		if(mers[id].hash >= minmers[insert_id].hash) 
			//		{
			//			break_id = insert_id + 1;
			//			break;
			//		}	
			//	}
			//	memcpy(minmers + break_id + 1, minmers + break_id, (l - break_id - 1) * sizeof(mer_info) );
			//	minmers[break_id] = mers[id];
			//}
			//if(pqueue.empty())
			//	pqueue.push(mers[id]);
			//else 
			if( pqueue.size() < l || mers[id].hash < pqueue.top().hash)
			{
				pqueue.push(mers[id]);
				if(pqueue.size() > l) pqueue.pop();
			}

		}

		lmers.clear();
		while(!pqueue.empty())
		{
			lmers.push_back(pqueue.top());
			pqueue.pop();
		}
		//std::sort(lmers.begin(), lmers.end(), [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });
		std::sort(minmers, minmers + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });
		//assert(lmers.size() == l);

		//	block(i, j, lmers[j].pos);
		for(unsigned j = 0; j < l; ++j)
		{
			memcpy(ptr, &seq.data()[minmers[j].pos], k);
			ptr += k;
		}
	}
	t2 = get_sec();
	std::cout << "omh main sketch time: " << t2 - t1 << std::endl;
}

double OrderMinHash::compare_sketches(const OSketch& sk1, const OSketch& sk2, ssize_t m, bool circular) {
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

double OrderMinHash::compare_sketch_pair(const char* p1, const char* p2, unsigned m, unsigned k, unsigned l, bool circular) {
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

double OrderMinHash::similarity(OrderMinHash & omh2){
	return compare_sketches(this->sk, omh2.getSektch());	
}

inline uint64_t hash_to_uint(const char * kmer, int k)
{
	uint8_t mask = 0x06; //FIXME: not general only works for DNA sequences, it's just a trick.
	uint64_t res = 0;
	for(int i = 0; i < k; i++)
	{
		uint8_t meri = (uint8_t)kmer[i];
		meri &= mask;
		meri >>= 1;
		res |= (uint64_t)meri;
		res <<= 2;
	}

	return res;
}	

}// namespace Sketch
