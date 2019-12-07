#ifndef _KMERSPECTRUM_
#define _KMERSPECTRUM_

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <vector>

using namespace std;

struct Bin{
	uint32_t BinID;
	double Frequency;
};

//int32_t JumpConsistentHash(uint64_t key, int32_t num_buckets);

void kmerSpectrumAddHash(vector <uint64_t> minimizerSketch, double * binsArr, int numBins);
void kmerSpectrumDump(double *binsArr, int numBins, vector <Bin> &kmerSpectrum);

	
//wrong logical implement: the minimizerSketch must be reflush;
//void getKmerSpectrum(vector <uint64_t> minimizerSketch, double * binsArr, int numBins, vector <Bin> &kmerSpectrum);



#endif








