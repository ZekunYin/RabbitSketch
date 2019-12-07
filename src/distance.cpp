#include "distance.h"

double getWJD(uint32_t *setA, uint32_t *setB, double *weightA, double *weightB, int sketchSize){
	double intersectElement = 0.0;
	double unionElement = 0.0;
	double distance = 1.0;

	for(int i = 0; i < sketchSize; i++){
		double curWeightA = abs(weightA[i]);
		double curWeightB = abs(weightB[i]);

		//get the intersection and union values
		if(setA[i] == setB[i]){
			if(curWeightA < curWeightB){
				intersectElement += curWeightA;
				unionElement += curWeightB;
			}
			else{
				intersectElement += curWeightB;
				unionElement += curWeightA;
			}
		}
		else{
			if(curWeightA > curWeightB){
				unionElement += curWeightA;
			}
			else{
				unionElement += curWeightB;
			}
		}

	}

	distance = 1 - intersectElement / unionElement;

	return distance;

}

				

		




