#ifndef _HISTOSKETCH_
#define _HISTOSKETCH_

#include <iostream>
#include <cstdio>
#include <stdint.h>

using namespace std;

double histoSketch_getSample(int i, int j, double freq, double * r, double * c, double * b, int sketchSize, int dimension);

//void histoSketchAddElement(uint64_t bin, double value, double * countMinSketch, int histoSketchLength, bool applyConceptDrift, double * r, double * c, double * b, int sketchSize, int dimension, uint32_t * histoSketch_sketch, double * histoSketch_sketchWeight);


void histoSketchAddElement(uint64_t bin, double value, double * countMinSketch, int histoSketchLength, bool applyConceptDrift, double decayWeight, double * r, double * c, double * b, int sketchSize, int dimension, uint32_t * histoSketch_sketch, double * histoSketch_sketchWeight);




#endif
