#include <stdio.h>
#include <iostream>
#include <string.h>
#include <list>
#include <cstring>
#include "Sketch.h"
#include "kseq.h"
#include "MurmurHash3.h"
#include "kseq.h"
#include "hash.h"

using namespace std;
using namespace Sketch;

//Writen by qzh.What I should do in MinHash constructor.
//Sketch::MinHash::MinHash(const Sketch::Parameters & parameters)
MinHash ::MinHash(Parameters parametersNew):parameters(parametersNew)
{
	minHashHeap = new MinHashHeap(parameters.use64, parameters.minHashesPerWindow, parameters.reads ?  parameters.minCov : 1);
	

}

void MinHash::update(char * seq)
{
	//addMinHashes(minHashHeap, seq, LENGTH, parameters);
	const uint64_t LENGTH = strlen(seq);
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
		for ( uint64_t i = 0; i < LENGTH; i++ )
		{
		    char base = seq[i];
		    
		    switch ( base )
		    {
		        case 'A': base = 'T'; break;
		        case 'C': base = 'G'; break;
		        case 'G': base = 'C'; break;
		        case 'T': base = 'A'; break;
		        default: break;
		    }
		    
		    seqRev[LENGTH - i - 1] = base;
		}
    }
    
    for ( uint64_t i = 0; i < LENGTH - kmerSize + 1; i++ )
    {
		// repeatedly skip kmers with bad characters
		bool bad = false;
		
		//Modified by qzh.To detect the correct of alphabet, but it consumes too much time. So we should optimize this process.
		//for ( uint64_t j = i; j < i + kmerSize && i + kmerSize <= LENGTH; j++ )
		//{
		//	if ( ! parameters.alphabet[seq[j]] )
		//	{
		//		i = j; // skip to past the bad character
		//		bad = true;
		//		break;
		//	}
		//}
		//
		if ( bad )
		{
			continue;
		}
		//	
		if ( i + kmerSize > LENGTH )
		{
			// skipped to end
			break;
		}
            
        const char *kmer_fwd = seq + i;
        const char *kmer_rev = seqRev + LENGTH - i - kmerSize;
        const char * kmer = (noncanonical || memcmp(kmer_fwd, kmer_rev, kmerSize) <= 0) ? kmer_fwd : kmer_rev;
        bool filter = false;
        
        hash_u hash = getHash(kmer, kmerSize, parameters.seed, parameters.use64);
        
		minHashHeap -> tryInsert(hash);
    }
    
    if ( ! noncanonical )
    {
        delete [] seqRev;
    }
}

void MinHash::printHashList()
{
	//setMinHashesForReference(reference, minHashHeap);
	HashList & hashlist = reference.hashesSorted;
	hashlist.clear();
	minHashHeap -> toHashList(hashlist);
	minHashHeap -> toCounts(reference.counts);
	hashlist.sort();

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
