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
#include <math.h>
#include <gsl/gsl_cdf.h>

//#ifdef USE_BOOST
//    #include <boost/math/distributions/binomial.hpp>
//    using namespace::boost::math;
//#else
//    #include <gsl/gsl_cdf.h>
//#endif



using namespace std;
using namespace Sketch;

//Writen by qzh.What I should do in MinHash constructor.
//Sketch::MinHash::MinHash(const Sketch::Parameters & parameters)
MinHash ::MinHash(Parameters parametersNew):parameters(parametersNew)
{
	minHashHeap = new MinHashHeap(parameters.use64, parameters.minHashesPerWindow, parameters.reads ?  parameters.minCov : 1);
	
	this->kmerSpace = pow(parameters.alphabetSize, parameters.kmerSize);
	//cerr << "kmerSpace init from pow is " << this->kmerSpace << endl;
	this->length = 0;

}

void MinHash::update(char * seq)
{
	//addMinHashes(minHashHeap, seq, LENGTH, parameters);
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

//	for(int i = 0; i < reference.hashesSorted.size(); i++){
//		if(parameters.use64)
//			cerr << "hash64 " <<  i << " " << reference.hashesSorted.at(i).hash64 << endl;
//		else
//			cerr << "hash32 " <<  i << " " << reference.hashesSorted.at(i).hash32 << endl;
//	}
//	return;
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
//	cout << "the size of hashesSortedRef is: " << hashesSortedRef.size() << endl;
//	cout << "the size of hashesSortedQry is: " << hashesSortedQry.size() << endl;
	cout << "the size of hashesSortedRef is: " << this->reference.hashesSorted.size() << endl;
	cout << "the size of hashesSortedQry is: " << msh->reference.hashesSorted.size() << endl;
	
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









