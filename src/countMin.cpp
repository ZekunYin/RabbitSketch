#include "countMin.h"
#include "jumpHash.h"

double countMinAdd(uint64_t element, double increment, bool applyScaling, double * countMinSketch){
	
	int g = ceil(2 / EPSILON);//math.h
	int d = ceil(log(1 - DELTA) / log(0.5));//math.h

	if(applyScaling){
		//TODO
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
	
	








