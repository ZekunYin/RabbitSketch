#ifndef __ORDERMINHASH_H__
#define __ORDERMINHASH_H__

#include <string>
#include <stdint.h>

namespace Sketch{

template<typename BT>
static void omh_pos(const std::string& seq, unsigned k, unsigned l, unsigned m, uint64_t mtSeed, BT block);

struct mer_info {
	size_t pos;
	uint64_t hash;
	unsigned occ;
	mer_info(size_t p, unsigned o, uint64_t h)
		: pos(p)
		  , hash(h)
		  , occ(o)
	{ }
};

inline uint64_t hash_to_uint(const char * kmer, int k);

}
#endif //__ORDERMINHASH_H__
