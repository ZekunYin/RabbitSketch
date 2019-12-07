#ifndef _COUNTMIN_
#define _COUNTMIN_
#include <iostream>
#include <stdint.h>
//#include <cstdio>
#include <cstdlib>
//#include <vector>
#include <cmath> //get the maxFloat64

using namespace std;

const double EPSILON = 0.001;
const double DELTA = 0.99;

//double countMinAdd(uint64_t element, double increment, bool applyScaling, double * countMinSketch, int d, int g);
double countMinAdd(uint64_t element, double increment, bool applyScaling, double * countMinSketch);











#endif
