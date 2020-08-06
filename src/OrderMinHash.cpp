
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
	omh_pos(seqStr, m_k, m_l, m_m, mtSeed, ptr);

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
	int lanes = 8;// vector lanes
	uint64_t buffer[lanes];

	for(unsigned i = 0; i < m; ++i) {
		const auto seed = gen64();//prg();
		size_t pend_size = (mers.size() / lanes) * lanes;
		//body
		for( int id = 0; id < pend_size; id+=lanes)
		{
			//uint64_t kmer_int = mers[id].int_hash;
			#pragma ivdep
			#pragma unroll(4)
			for(int vid = 0; vid < lanes; vid++)	
			{
				uint64_t kmer_int = intHash[id + vid];
		 	    kmer_int += occ[id + vid] * weight;
				buffer[vid] ^= seed; 
    			buffer[vid] ^= buffer[vid] >> 33;
    			buffer[vid] *= 0xff51afd7ed558ccd;
    			buffer[vid] ^= buffer[vid] >> 33;
    			buffer[vid] *= 0xc4ceb9fe1a85ec53;
    			buffer[vid] ^= buffer[vid] >> 33;

			}
			//meri.hash = hash.digest();
			for(int vid = 0; vid < lanes; vid++)	
			{
				mers[id + vid].hash = buffer[vid];
				if(pqueue.empty())
					pqueue.push(mers[id + vid]);
				else if(mers[id + vid].hash < pqueue.top().hash || pqueue.size() < l)
				{
					pqueue.push(mers[id + vid]);
					if(pqueue.size() > l) pqueue.pop();
				}
			}
		}

		//tail
		for( int id = pend_size; id < pend_size; id++)
		{

			uint64_t kmer_int = intHash[id];
		    kmer_int += occ[id] * weight;
			mers[id].hash = mc::murmur3_fmix(kmer_int, seed);
			if(pqueue.empty())
				pqueue.push(mers[id]);
			else if(mers[id].hash < pqueue.top().hash || pqueue.size() < l)
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
		std::sort(lmers.begin(), lmers.end(), [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });
		assert(lmers.size() == l);

		//	block(i, j, lmers[j].pos);
		for(unsigned j = 0; j < l; ++j)
		{
			memcpy(ptr, &seq.data()[lmers[j].pos], k);
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
