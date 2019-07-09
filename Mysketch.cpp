#include <stdio.h>
#include <iostream>
#include <string.h>
#include <list>
#include "Sketch.h"
#include "kseq.h"
#include "MurmurHash3.h"
#include "kseq.h"
#include "hash.h"

using namespace std;


void sketchOutput(Sketch::SketchInput * input)
{
	const Sketch::Parameters & parameters = input->parameters;
	
	Sketch::SketchOutput * output = new Sketch::SketchOutput();

    cerr << "The seq is : " << input->seq << endl;
    cerr << "The length is : " << input->length << endl;

	//parameters.noncanonical = false;
    //parameters.preserveCase = false;    
    //parameters.use64 = false;

	output->references.resize(1);
	Sketch::Reference & reference = output->references[0];
	
	//Modified by qzh. I have remove the name and comment from the SketchInput.
	//reference.name = input->name;
	//reference.comment = input->comment;
	reference.length = input->length;
	reference.hashesSorted.setUse64(parameters.use64);
	
	MinHashHeap minHashHeap(parameters.use64, parameters.minHashesPerWindow, parameters.reads ? parameters.minCov : 1);
    addMinHashes(minHashHeap, input->seq, input->length, parameters);
	setMinHashesForReference(reference, minHashHeap);

	for(int i = 0; i < reference.hashesSorted.size(); i++){
		if(parameters.use64)
			cerr << "hash64 " <<  i << " " << reference.hashesSorted.at(i).hash64 << endl;
		else
			cerr << "hash32 " <<  i << " " << reference.hashesSorted.at(i).hash32 << endl;
	}

	//return output;
}

//addMinHashes
void addMinHashes(MinHashHeap & minHashHeap, char * seq, uint64_t length, const Sketch::Parameters & parameters)
{
    int kmerSize = parameters.kmerSize;
	//cerr << "kmerSize is = "<< kmerSize << endl;
    uint64_t mins = parameters.minHashesPerWindow;
	//cerr << "mins is = "<< mins << endl;
    bool noncanonical = parameters.noncanonical;//False.
	//cerr << "noncanonical is = "<< noncanonical << endl;
    
    // Determine the 'mins' smallest hashes, including those already provided
    // (potentially replacing them). This allows min-hash sets across multiple
    // sequences to be determined.
    
    // uppercase TODO: alphabets?
    //
    for ( uint64_t i = 0; i < length; i++ )
    {
        if ( ! parameters.preserveCase && seq[i] > 96 && seq[i] < 123 )
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
    
    for ( uint64_t i = 0; i < length - kmerSize + 1; i++ )
    {
		// repeatedly skip kmers with bad characters
		//
		bool bad = false;
		
		//Modified by qzh.To detect the correct of alphabet, but it consumes too much time. So we should optimize this process.
		//for ( uint64_t j = i; j < i + kmerSize && i + kmerSize <= length; j++ )
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
		if ( i + kmerSize > length )
		{
			// skipped to end
			break;
		}
            
        const char *kmer_fwd = seq + i;
        const char *kmer_rev = seqRev + length - i - kmerSize;
        const char * kmer = (noncanonical || memcmp(kmer_fwd, kmer_rev, kmerSize) <= 0) ? kmer_fwd : kmer_rev;
        bool filter = false;
        
        hash_u hash = getHash(kmer, kmerSize, parameters.seed, parameters.use64);
        
		minHashHeap.tryInsert(hash);
    }
    
    if ( ! noncanonical )
    {
        delete [] seqRev;
    }
}

void reverseComplement(const char * src, char * dest, int length)
{
    for ( int i = 0; i < length; i++ )
    {
        char base = src[i];
        
        switch ( base )
        {
            case 'A': base = 'T'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
            case 'T': base = 'A'; break;
            default: break;
        }
        
        dest[length - i - 1] = base;
    }
}

void setMinHashesForReference(Sketch::Reference & reference, const MinHashHeap & hashes)
{
    HashList & hashList = reference.hashesSorted;
    hashList.clear();
    hashes.toHashList(hashList);
    hashes.toCounts(reference.counts);
    hashList.sort();
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

