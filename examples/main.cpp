#include "Sketch.h"
#include <iostream>
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"
#include <vector>
#include <math.h>
#include <random>

using namespace std;

KSEQ_INIT(gzFile, gzread)

double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

//static void getCWS(double *r, double *c, double *b, int sketchSize, int dimension){
////	cerr << "successful malloc r, c, b in getCWS" << endl;
//	const int DISTRIBUTION_SEED = 1;
//    default_random_engine generator(DISTRIBUTION_SEED);
//    gamma_distribution<double> gamma(2.0,1.0);
//    uniform_real_distribution<double> uniform(0.0,1.0);
//
//    for (int i = 0; i < sketchSize * dimension; ++i){
//        r[i] = gamma(generator);
//        c[i] = log(gamma(generator));
//        b[i] = uniform(generator) * r[i];
//    }
//}	


int main(int argc, char* argv[])
{

	gzFile fp1;
	kseq_t *ks1;
	
	fp1 = gzopen(argv[1],"r");
	if(NULL == fp1){
		fprintf(stderr,"Fail to open file: %s\n", argv[1]);
		return 0;
	}


	vector<string> vseq;
	vector<int> vlength;
	int count = 0;

	Sketch::WMHParameters parameter;
	parameter.kmerSize = 21;
	parameter.sketchSize = 50;
	parameter.windowSize = 20;
	parameter.r = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
	parameter.c = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
	parameter.b = (double *)malloc(parameter.sketchSize * pow(parameter.kmerSize, 4) * sizeof(double));
	getCWS(parameter.r, parameter.c, parameter.b, parameter.sketchSize, pow(parameter.kmerSize, 4));
	vector<Sketch::WMinHash *> vwmh; 
	vector<Sketch::MinHash *> vmh; 
	vector<Sketch::OrderMinHash > vomh; 
	vector<Sketch::HyperLogLog> vhlog;
	ks1 = kseq_init(fp1);
	while(1){
		int length = kseq_read(ks1);
		if(length < 0){
			break;
		}

		Sketch::WMinHash * wmh1 = new Sketch::WMinHash(parameter);
		Sketch::MinHash * mh1 = new Sketch::MinHash();
		Sketch::OrderMinHash omh1;
    	static const size_t BITS = 20; //24
		Sketch::HyperLogLog t(BITS);

		cerr << "end the wmh construction" << endl;
		wmh1->update(ks1->seq.s);
		mh1->update(ks1->seq.s);	
		omh1.buildSketch(ks1->seq.s);
		t.update(ks1->seq.s);
		cerr << "end the wmh update" << endl;
		vwmh.push_back(wmh1);
		vmh.push_back(mh1);
		vomh.push_back(omh1);
		vhlog.push_back(t);
		cerr << "the index of the wmh is: " << count << endl;
		count++;

	}


//	for(int i = 0; i < count; i++){
//		//Sketch::WMinHash wmh = new Sketch::WMinHash();
//	}
	cout << "begin to compute the WMH distance: "  << endl;
	printf("=====================================\t WMinHash \t MinHash \t OMinHash \t HyperLog\n");
	for(int i = 0; i < count; i++){
		for(int j = i+1; j < count; j++){
			double distance0 = vwmh[i]->distance(vwmh[j]);
			double distance1 = 1.0 - vmh[i]->jaccard(vmh[j]);
			double distance2 = vomh[i].distance(vomh[j]);
			double distance3 = vhlog[i].distance(vhlog[j]);
			printf("the distance of seq[%d] and seq[%d] is:\t %lf \t %lf \t %lf \t %lf\n", i, j, distance0, distance1, distance2, distance3);
		}
	}
	cout << "end to compute the WMH distance;"  << endl;



	return 0;
	
}

