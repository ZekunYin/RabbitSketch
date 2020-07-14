// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef hash_h
#define hash_h

#include <inttypes.h>
#include "MurmurHash3.h"

typedef uint32_t hash32_t;
typedef uint64_t hash64_t;

union hash_u
{
    hash32_t hash32;
    hash64_t hash64;
};

enum hash_type {
    type_MurmurHash3_x86_32, type_MurmurHash3_x86_128, type_MurmurHash3_x64_128  
};

void (*hash_array[])(const void *, int, uint32_t, void *) = 
    {MurmurHash3_x86_32, MurmurHash3_x86_128, MurmurHash3_x64_128};

hash_u getHash(const char * seq, int length, uint32_t seed, bool use64);
bool hashLessThan(hash_u hash1, hash_u hash2, bool use64);

void sketchHashAlgo(const void * key, int len, uint32_t seed, void * out,
                    void (*pHashAlgo)(const void *, int, uint32_t, void *))
    {pHashAlgo(key, len, seed, out);};

#endif
