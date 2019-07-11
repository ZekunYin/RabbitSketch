// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

//Modified by qzh.
#ifndef Sketch_h
#define Sketch_h

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include "MinHashHeap.h"
#include "ThreadPool.h"

namespace Sketch{

typedef uint64_t hash_t;

struct Parameters
{
	Parameters()
		:
			parallelism(1),
			kmerSize(0),
			alphabetSize(0),
			preserveCase(false),
			use64(false),
			seed(0),
			error(0),
			warning(0),
			minHashesPerWindow(0),
			windowSize(0),
			windowed(false),
			concatenated(false),
			noncanonical(false),
			reads(false),
			memoryBound(0),
			minCov(1),
			targetCov(0),
			genomeSize(0)
	{
		memset(alphabet, 0, 256);
	}

	Parameters(const Parameters & other)
		:
			parallelism(other.parallelism),
			kmerSize(other.kmerSize),
			alphabetSize(other.alphabetSize),
			preserveCase(other.preserveCase),
			use64(other.use64),
			seed(other.seed),
			error(other.error),
			warning(other.warning),
			minHashesPerWindow(other.minHashesPerWindow),
			windowSize(other.windowSize),
			windowed(other.windowed),
			concatenated(other.concatenated),
			noncanonical(other.noncanonical),
			reads(other.reads),
			memoryBound(other.memoryBound),
			minCov(other.minCov),
			targetCov(other.targetCov),
			genomeSize(other.genomeSize)
	{
		memcpy(alphabet, other.alphabet, 256);
	}

	int parallelism;
	int kmerSize;
	bool alphabet[256];
	uint32_t alphabetSize;
	bool preserveCase;
	bool use64;
	uint32_t seed;
	double error;
	double warning;
	uint64_t minHashesPerWindow;
	uint64_t windowSize;
	bool windowed;
	bool concatenated;
	bool noncanonical;
	bool reads;
	uint64_t memoryBound;
	uint32_t minCov;
	double targetCov;
	uint64_t genomeSize;
};

struct PositionHash
{
	PositionHash(uint32_t positionNew, hash_t hashNew) :
		position(positionNew),
		hash(hashNew)
	{}

	uint32_t position;
	hash_t hash;
};

struct Locus
{
	Locus(uint32_t sequenceNew, uint32_t positionNew)
		:
			sequence(sequenceNew),
			position(positionNew)
	{}

	uint32_t sequence;
	uint32_t position;
};

typedef std::unordered_set<hash_t> Hash_set;

struct Reference
{
	// no sequence for now
	std::string name;
	std::string comment;
	uint64_t length;
	HashList hashesSorted;
	std::vector<uint32_t> counts;
};

class MinHash
{

	public:
		//Modified by qzh.Remove const and the length in the function.	
		//MinHash(const Parameters & parameters);
		//void addMinHashes(char * seq, uint64_t length);
		//void update(char * seq, uint64_t length);
		MinHash(Parameters parametersNew);
		void update(char * seq);

		double jaccard(MinHash & msh);			
		double dist(MinHash & msh);

		HashList & getHashList();
		void printHashList();
		void writeToFile();

	private:
		//Modified by qzh.Define the seq.
		char * seq;
		double kmerSpace;

		Parameters parameters;
		MinHashHeap * minHashHeap;
		Reference reference;

};

class WMinHash{
};

class OMinHash{
};
}//namespace sketch

#endif //Sketch_h
