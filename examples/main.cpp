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

void getCWS(double *r, double *c, double *b, int sketchSize, int dimension){
//	cerr << "successful malloc r, c, b in getCWS" << endl;
	const int DISTRIBUTION_SEED = 1;
    default_random_engine generator(DISTRIBUTION_SEED);
    gamma_distribution<double> gamma(2.0,1.0);
    uniform_real_distribution<double> uniform(0.0,1.0);

    for (int i = 0; i < sketchSize * dimension; ++i){
        r[i] = gamma(generator);
        c[i] = log(gamma(generator));
        b[i] = uniform(generator) * r[i];
    }
}	


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


//	double distance;
//	double time1 = get_sec();	

//	//Sketch::WMinHash * wmh1 = new Sketch::WMinHash(21, 50, 9, 0.0);
//	//Sketch::WMinHash * wmh2 = new Sketch::WMinHash(21, 50, 9, 0.0);
//	Sketch::WMinHash * wmh1 = new Sketch::WMinHash();
//	Sketch::WMinHash * wmh2 = new Sketch::WMinHash();
//	//wmh1->setHistoSketchSize(500);
//	//wmh2->setHistoSketchSize(500);
//	wmh1->update(seq1);
//	wmh2->update(seq2);
//
//	double time3 = get_sec();
//
//	distance = wmh1->distance(wmh2);
//
//	double time2 = get_sec();
//	cout << "Algorithm\t" << "distance\t" << "totaltime\t" << "sketchtime\t" << endl;
//	cout << "WMinHash\t"  << distance << "\t" << time2 - time1 << "\t" << time3 - time1 << endl;
//	//cout << "WMinHash: the distance(1-WJ) is: " << distance << endl;
//	//cout << "WMinHash time: " << time2 - time1  << endl;
//	//cout << "WMinHash sketch time: " << time3 - time1 << endl;
//	//cout << "=======================================================" << endl;
//	
//	time1 = get_sec();
//
//	Sketch::MinHash * mh1 = new Sketch::MinHash();
//	Sketch::MinHash * mh2 = new Sketch::MinHash();
//	mh1->update(seq1);	
//	//std::cout << "mh1 set size: " << mh1->count() << std::endl;
//	mh2->update(seq2);	
//	//std::cout << "mh2 set size: " << mh2->count() << std::endl;
//
//	time3 = get_sec();
//
//	double minhash_jac = mh1->jaccard(mh2);	
//	//double minhash_jac = mh1->mdistance(mh2);	
//
//	time2 = get_sec();
//
//	cout << "MinHash\t"  << 1.0 - minhash_jac << "\t" << time2 - time1 << "\t" << time3 - time1 << endl;
//	//cout << "MinHash: the distance(1-J) is: " << 1.0 - minhash_dist << endl;
//	//cout << "MinHash time: " << time2 - time1  << endl;
//	//cout << "MinHash sketch time: " << time3 - time1  << endl;
//	//cout << "=======================================================" << endl;
//
//
//	time1  = get_sec();
//
//	Sketch::OrderMinHash omh1;
//	Sketch::OrderMinHash omh2;
//	omh1.buildSketch(seq1);
//	omh2.buildSketch(seq2);
//
//	time3 = get_sec();
//
//	double odist = omh1.distance(omh2);
//
//	time2 = get_sec();
//
//	cout << "OMinHash\t"  << odist  << "\t" << time2 - time1 << "\t" << time3 - time1 << endl;
//	//cout << "OMinHash: the distance(1-J) is: " << odist << endl;
//	//cout << "OMinHash time: " << time2 -  time1 << endl;
//	//cout << "OMinHash sketch: " << time3 -  time1 << endl;
//	//cout << "=======================================================" << endl;
//
//
//	time1  = get_sec();
//
//    static const size_t BITS = 20; //24
//	Sketch::HyperLogLog t(BITS);
//	Sketch::HyperLogLog t1(BITS);
//
//    t.update(seq1);
//    t1.update(seq2);
//
//	time3 = get_sec();
//
//    double dist = t1.distance(t);
//
//	time2  = get_sec();
//
//	cout << "HLL\t"  << 1.0 - dist  << "\t" << time2 - time1 << "\t" << time3 - time1 << endl;
//    //cout << "HLL distance(1-J) = " <<  1.0 - dist << endl;
//	//cout << "HLL time: " << time2 -  time1 << endl;
//	//cout << "HLL sketch time: " << time3 -  time1 << endl;
//
//	kseq_destroy(ks1);
//	gzclose(fp1);

	return 0;
	
}

