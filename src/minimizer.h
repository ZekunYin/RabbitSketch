#ifndef _MINIMIZER_
#define _MINIMIZER_

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <stdint.h>

using namespace std;
struct Pair{
	uint64_t X;//the hash value.
	int Y;//the location of the minimizer in the sequence.
};

//vector <uint64_t> findMinimizers(int k, int w, string s);
void findMinimizers(int k, int w, string s, vector <uint64_t> &minimizerSketch);

uint64_t hash64(uint64_t key, uint64_t mask);
bool findElement(vector <uint64_t> arr, uint64_t element);





#endif
