#include "histoSketch.h"
#include <math.h>


int32_t JumpConsistentHash(uint64_t key, int32_t num_buckets){
	int64_t b = -1, j = 0;
	while (j < num_buckets) {
		b = j;
		key = key * 2862933555777941757ULL + 1;
		j = (b + 1) * (double(1LL << 31) / double((key >> 33) + 1));
	}
	return b;
}

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

void countMinScale(double * countMinSketch, double decayWeight){

	int g = ceil(2 / EPSILON);
	int d = ceil(log(1 - DELTA) / log(0.5));

	for(int i = 0; i < d; i++){
		for(int j = 0; j < g; j ++){
			countMinSketch[i * g + j] = countMinSketch[i * g + j] * decayWeight;
		}
	}

}

double countMinAdd(uint64_t element, double increment, bool applyScaling, double decayWeight, double * countMinSketch){
	
	int g = ceil(2 / EPSILON);//math.h
	int d = ceil(log(1 - DELTA) / log(0.5));//math.h

	if(applyScaling){
		//TODO
		//cout << "do the countMinScale" << endl;
		countMinScale(countMinSketch, decayWeight);
		//countMinSketch.scale();
	}

	//double currentMinimum = HUGE_VALL;//use c++11: huge long double value
	double currentMinimum = 10000.0;//use c++11: huge long double value
	//cout << currentMinimum << endl;
	
	for(int i = 0; i < d; i++){
		uint64_t hash = element + (uint64_t)i * element;

		int j = JumpConsistentHash(hash, g);
	//	cout << j << endl;

		if(increment != 0.0){//actrually don't need the if
			countMinSketch[i * g + j] += increment;
		}

		if(countMinSketch[i * g + j] < currentMinimum){
			currentMinimum = countMinSketch[i * g + j];
		}
	
	}

	return currentMinimum;
}
	

//need to be implement.//TODO the dimention of r, c, b.
double histoSketch_getSample(int i, int j, double freq, double * r, double * c, double * b, int sketchSize, int dimension){
	//the b[j][i] need to be relocated the index.
	if(j >= sketchSize){
		cerr <<"out of bound sketchSize!" << endl;
		exit(1);
	}
	if(i >= dimension){
		cerr << "out of bound dimension!" << endl;
		exit(1);
	}

	//double yka = exp(log(freq) - b[j][i]);
//	printf("the b is: %lf\n", b[j * dimension + i]);
	double yka = exp(log(freq) - b[j * dimension + i]);
	//double result = c[j][i] / (yka * exp(r[j][i]));
	double result = c[j * dimension + i] / (yka * exp(r[j * dimension + i]));
	return result;
}



void histoSketchAddElement(uint64_t bin, double value, double * countMinSketch, int histoSketchLength, bool applyConceptDrift, double decayWeight, double * r, double * c, double * b, int sketchSize, int dimension, uint32_t * histoSketch_sketch, double * histoSketch_sketchWeight){
	
	double estiFreq = countMinAdd(bin, value, applyConceptDrift, decayWeight, countMinSketch);

	for(int i = 0; i < histoSketchLength; i++){
		//get the CWS value(A_Ka) for the incoming element
		double aka = histoSketch_getSample((int)bin, i, estiFreq, r, c, b, sketchSize, dimension);	
//		printf("the aka is: %lf\n", aka);
//		exit(1);

		//get the current minimum in the histosketchSlot, accounting for concept drift if requrested
		double curMin;
		if(applyConceptDrift){
			//TODO
			curMin = histoSketch_sketchWeight[i] /decayWeight;
		}
		else{
			curMin = histoSketch_sketchWeight[i];
		}

		if(aka < curMin){
			histoSketch_sketch[i] = (uint32_t)bin;
			histoSketch_sketchWeight[i] = aka;
		}
	}

}
		












		
