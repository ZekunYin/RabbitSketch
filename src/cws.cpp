#include "cws.h"
#include <random>
#include <cmath>
//#include "readFile.h" //for tmp: read file for get r, c, b;

using namespace std;

//get the cws by reading file.
//void getCWS(string fileName, double *rcb, int sketchSize, int dimension){
//	ifstream data(fileName.c_str());
//	string line;
//	if(!data.is_open()){
//		cerr << "unable open the file!" << endl;
//		exit(1);
//	}
//	long curIndex = 0;
//
//	while(!data.eof()){
//		getline(data, line);
//		stringstream s1;
//		s1 << line;
//		double theNumber;
//		s1 >> theNumber;
//
//		if(curIndex < sketchSize * dimension){
//
//			rcb[curIndex] = theNumber;
//			curIndex++;
//		}
//	}
//
//}

const int DISTRIBUTION_SEED = 1;

void getCWS(double *r, double *c, double *b, int sketchSize, int dimension){
    default_random_engine generator(DISTRIBUTION_SEED);
    gamma_distribution<double> gamma(2.0,1.0);
    uniform_real_distribution<double> uniform(0.0,1.0);

    for (int i = 0; i < sketchSize * dimension; ++i){
        r[i] = gamma(generator);
        c[i] = log(gamma(generator));
        b[i] = uniform(generator) * r[i];
    }
}	






