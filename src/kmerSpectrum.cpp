//#include <iostream>
//#include <cstdio>
//#include <cstdlib>
//#include <stdint.h>
//#include <vector>
#include "kmerSpectrum.h"
#include "jumpHash.h"

using namespace std;

//struct Bin{
//	uint32_t BinID;
//	double Frequency;
//};
	
//int32_t JumpConsistentHash(uint64_t key, int32_t num_buckets){
//	int64_t b = -1, j = 0;
//	while (j < num_buckets) {
//		b = j;
//		key = key * 2862933555777941757ULL + 1;
//		j = (b + 1) * (double(1LL << 31) / double((key >> 33) + 1));
//	}
//	return b;
//}


void kmerSpectrumAddHash(vector <uint64_t> minimizerSketch, double * binsArr, int numBins){
	for(int i = 0; i < minimizerSketch.size(); i++){
		
		uint32_t bin = JumpConsistentHash(minimizerSketch[i], numBins);
//		cout << "the bin is: " << bin << endl;

		//record the bin id; check the kmerSpectrum.bv.add ? 
		{

		}

		binsArr[bin] += 1.0;

	}
}

void kmerSpectrumDump(double *binsArr, int numBins, vector <Bin> &kmerSpectrum){
	//implement the KmerSpectrum.Dump in go.
	for(int i = 0; i < numBins; i++){
		if(binsArr[i] != 0.0){
			Bin binDemo;
			binDemo.BinID = i;
			binDemo.Frequency = binsArr[i];
			kmerSpectrum.push_back(binDemo);
		}
	}

}
	
