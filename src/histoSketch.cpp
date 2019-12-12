#include "histoSketch.h"
#include "countMin.h"



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
		












		
