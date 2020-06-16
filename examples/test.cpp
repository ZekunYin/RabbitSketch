#include "Sketch.h"
#include <iostream>
#include <sys/time.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
using namespace std;
using namespace Sketch;

double get_sec(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}


int main(int argc, char* argv[])
{
	int k = 21;
	Sketch::Parameters parameters;
	
    parameters.alphabetSize = 4;
    parameters.kmerSize = k;
	parameters.use64 = pow(parameters.alphabetSize, parameters.kmerSize) > pow(2, 32);
	parameters.seed = 42;
	parameters.minHashesPerWindow  = 1000;


	parameters.set_numBins(pow(parameters.kmerSize, parameters.alphabetSize));
	parameters.set_minimizerWindowSize(9);
	parameters.set_histoSketch_sketchSize(60);
	parameters.set_histoSketch_dimension(194481);
	parameters.set_paraDecayWeight(0.1);


//	parameters.numBins = pow(parameters.kmerSize, parameters.alphabetSize);
//	parameters.minimizerWindowSize = 9;
//	parameters.histoSketch_sketchSize = 500;
//	parameters.histoSketch_dimension = 194481;
//	parameters.paraDecayWeight = 0.1;

	gzFile fp1;
	gzFile fp2;
	kseq_t *ks1;
	kseq_t *ks2;
	
	fp1 = gzopen(argv[1],"r");
	fp2 = gzopen(argv[2],"r");
	if(NULL == fp1){
		fprintf(stderr,"Fail to open file: %s\n", argv[1]);
		return 0;
	}

	if(NULL == fp2){
		fprintf(stderr,"Fail to open file: %s\n", argv[2]);
		return 0;
	}

	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);
/*	
	int nums = 100;
	if(argc >=4)
		nums = stoi(argv[3]);
	int count = 0;

	std::vector<HyperLogLog*> vect_hll;
    static const size_t BITS = 10; //24
	
	//readfile & sketch
	while( kseq_read(ks) >= 0 ){
		if(argc >= 4 && count+1 > nums )
			break;
		
		HyperLogLog * test = new HyperLogLog(BITS);
		test->update(ks->seq.s);
		vect_hll.push_back(test);
	
		count++;
	}
	
	//distance
	for(int i=0; i<vect_hll.size(); ++i) {
		fprintf(stdout, "current is %d \n",i);
		for(int j=i; j<vect_hll.size(); ++j){
			double dist = vect_hll[i]->distance(*vect_hll[j]);
			double esti = (vect_hll[i]->merge(*vect_hll[j])).creport();
			fprintf(stdout,"Seq[%d] and seq[%d] distance(J): %lf , estimation = %lf \n", i, j, dist, esti);
			//fprintf(stdout,"Distance between sketch[%d] and sketch[%d]: %lf \n", i, j, dist);
		}
	}
*/
	int l1 = kseq_read(ks1);
	int l2 = kseq_read(ks2);
	char * seq1 = ks1->seq.s;
	char * seq2 = ks2->seq.s;

	double distance;
	double time1 = get_sec();	

	Sketch::WMinHash * wmh1 = new Sketch::WMinHash(parameters);
	Sketch::WMinHash * wmh2 = new Sketch::WMinHash(parameters);
	cerr << "start the update seq1" << endl;
	wmh1->update(seq1);
	cerr << "start the update seq2" << endl;
	wmh2->update(seq2);
	cerr << "before getWMinHash" << endl;
//	wmh1->getWMinHash();
//	wmh2->getWMinHash();
	cerr << "before distance " << endl;
	


	double time3 = get_sec();

	distance = wmh1->distance(wmh2);
	double time2 = get_sec();
	cout << "WMinHash: the distance(1-WJ) is: " << distance << endl;
	cout << "WMinHash time: " << time2 - time1  << endl;
	cout << "WMinHash sketch time: " << time3 - time1 << endl;
	
	time1 = get_sec();
	Sketch::MinHash * mh1 = new Sketch::MinHash(parameters);
	Sketch::MinHash * mh2 = new Sketch::MinHash(parameters);
	mh1->update(seq1);	
	mh2->update(seq2);	
	time3 = get_sec();

	double minhash_dist = mh1->jaccard(mh2);	

	time2 = get_sec();
	cout << "MinHash: the distance(1-J) is: " << 1.0 - minhash_dist << endl;
	cout << "MinHash time: " << time2 - time1  << endl;
	cout << "MinHash sketch time: " << time3 - time1  << endl;


	time1  = get_sec();
	Sketch::OMinHash omh1(parameters, seq1);
	Sketch::OMinHash omh2(parameters, seq2);
	time3 = get_sec();
	double odist = omh1.distance(omh2);
	time2 = get_sec();
	cout << "OMinHash: the distance(1-J) is: " << odist << endl;
	cout << "OMinHash time: " << time2 -  time1 << endl;
	cout << "OMinHash sketch: " << time3 -  time1 << endl;


	time1  = get_sec();
    static const size_t BITS = 20; //24
	HyperLogLog t(BITS);
	HyperLogLog t1(BITS);

    t.update(seq1);
    t1.update(seq2);
	time3 = get_sec();

    double dist = t1.distance(t);
	time2  = get_sec();
    cout << "HLL distance(1-J) = " <<  1.0 - dist << endl;
	cout << "HLL time: " << time2 -  time1 << endl;
	cout << "HLL sketch time: " << time3 -  time1 << endl;

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);

	return 0;
	
}

