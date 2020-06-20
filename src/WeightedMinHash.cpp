#include "Sketch.h"
#include <random>
#include "histoSketch.h"

namespace Sketch{

//start the WMinHash@xxm
WMinHash::WMinHash(Parameters parametersNew)//:parameters(parametersNew)
{	
	parameters = parametersNew;

	binsArr = (double *) malloc (parameters.get_numBins() * sizeof(double));
	for(int i = 0; i < parameters.get_numBins(); i++){
		binsArr[i] = 0.0;
	}

	int g = ceil(2 / EPSILON);
	int d = ceil(log(1 - DELTA) / log(0.5));
	countMinSketch = (double *) malloc (d * g *sizeof(double));
	for(int i = 0; i < d * g; i++){
		countMinSketch[i] = 0.0;
	}

	//r = (double *) malloc (parameters.get_histoSketch_sketchSize() * parameters.get_histoSketch_dimension() * sizeof(double));
	//c = (double *) malloc (parameters.get_histoSketch_sketchSize() * parameters.get_histoSketch_dimension() * sizeof(double));
	//b = (double *) malloc (parameters.get_histoSketch_sketchSize() * parameters.get_histoSketch_dimension() * sizeof(double));
	
	//getCWS(r, c, b, parameters.get_histoSketch_sketchSize(), parameters.get_histoSketch_dimension());
	//double * r = parameters.getR();
	//double * c = parameters.getC();
	//double * b = parameters.getB();
	//		printf("the point in WMinHash of r is: %p\n",r);
	//		printf("the point in WMinHash of c is: %p\n",c);
	//		printf("the point in WMinHash of b is: %p\n",b);
//	for(int i = 0; i < parameters.get_histoSketch_sketchSize() * parameters.get_histoSketch_dimension(); i++){
//		cout << "the r[" << i << "] is: " << parameters.getR()[i] << endl;
//		cout << "the c[" << i << "] is: " << parameters.getC()[i] << endl;
//		cout << "the b[" << i << "] is: " << parameters.getB()[i] << endl;
//	}
	cerr << "the paremeters.get_histoSketch_sketchSize() in WMinHash is: " << parameters.get_histoSketch_sketchSize() << endl;
	
	histoSketch_sketch = (uint32_t *) malloc (parameters.get_histoSketch_sketchSize() * sizeof(uint32_t));
	histoSketch_sketchWeight = (double *) malloc (parameters.get_histoSketch_sketchSize() * sizeof(double));
	
	//add the applyConceptDrift and decayWeight.
	if(parameters.get_paraDecayWeight() < 0.0 || parameters.get_paraDecayWeight() > 1.0){
		cerr << "the paraDecayWeight must between 0.0 and 1.0 " << endl;
		exit(1);
	}
	else{
		applyConceptDrift = true;
	}
	if(parameters.get_paraDecayWeight() == 1.0){
		applyConceptDrift = false;
	}
	decayWeight = 1.0;
	if(applyConceptDrift){
		decayWeight = exp(-parameters.get_paraDecayWeight());
	}


	needToCompute = true;

}

void WMinHash::update(char * seq)
{
	int k = parameters.kmerSize;
	int w = parameters.get_minimizerWindowSize();
	int numBins = parameters.get_numBins();
	int histoSketch_sketchSize = parameters.get_histoSketch_sketchSize();
	int histoSketch_dimension = parameters.get_histoSketch_dimension();
	cerr << "in the update the histoSketch_sketchSize is: " << parameters.get_histoSketch_sketchSize() << endl;

	findMinimizers(k, w, seq, sketches);
	kmerSpectrumAddHash(sketches, binsArr, parameters.get_numBins());
//	for(int i = 0; i < parameters.numBins; i++){
//		if(binsArr[i] != 0){
//			printf("%d\t%lf\n", i, binsArr[i]);
//		}
//	}
	sketches.clear();
	needToCompute = true;

}

void WMinHash::computeHistoSketch()
{
	kmerSpectrumDump(binsArr, parameters.get_numBins(), kmerSpectrums);
	//cerr << "finish the kmerSpectrumDump " << endl;
	int hsSketchSize = parameters.get_histoSketch_sketchSize();
	int hsDimension = parameters.get_histoSketch_dimension();
	//cerr << "the hsSketchSize is " << hsSketchSize << endl;
	//cerr << "the hsDimension is " << hsDimension << endl;

	for(int i = 0; i < kmerSpectrums.size(); i++){
	//cerr << "finish the " << i << " iterator to addElement" << endl;
		//histoSketchAddElement((uint64_t)kmerSpectrums[i].BinID, kmerSpectrums[i].Frequency, countMinSketch, parameters.get_histoSketch_sketchSize(), applyConceptDrift, decayWeight, r, c, b, parameters.get_histoSketch_sketchSize(), parameters.get_histoSketch_dimension(), histoSketch_sketch, histoSketch_sketchWeight);
		histoSketchAddElement((uint64_t)kmerSpectrums[i].BinID, kmerSpectrums[i].Frequency, countMinSketch, hsSketchSize, applyConceptDrift, decayWeight, parameters.getR(), parameters.getC(), parameters.getB(), hsSketchSize, hsDimension, histoSketch_sketch, histoSketch_sketchWeight);
	}
	//cerr << "finish the histoSketchAddElement " << endl;

}

void WMinHash::getWMinHash(){
	if(needToCompute){
		computeHistoSketch();
		needToCompute = false;
	}
	cout << "the sketch is: " << endl;
	for(int i = 0; i < parameters.get_histoSketch_sketchSize(); i++){
		printf("%d\n", histoSketch_sketch[i]);
	}

	cout << "the sketchweight is: " << endl;
	for(int i = 0; i < parameters.get_histoSketch_sketchSize(); i++){
		printf("%lf\n", histoSketch_sketchWeight[i]);
	}

}

double WMinHash::wJaccard(WMinHash * wmh)
{
	cout << "the needToCompute is: " << needToCompute << endl;
	if(needToCompute){
		computeHistoSketch();
		needToCompute = false;
	}
	if(wmh->needToCompute){
		wmh->computeHistoSketch();
		wmh->needToCompute = false;
	}
	double intersectElement = 0.0;
	double unionElement = 0.0;
	double jaccard = 0.0;
	cerr << "the parameters.get_histoSketch_sketchSize() is: " << parameters.get_histoSketch_sketchSize() << endl;
	for(int i = 0 ; i < parameters.get_histoSketch_sketchSize(); i++){

		double curWeightA = abs(histoSketch_sketchWeight[i]);
		double curWeightB = abs(wmh->histoSketch_sketchWeight[i]);

		//get the intersection and union values
		if(histoSketch_sketch[i] == wmh->histoSketch_sketch[i]){
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
	jaccard = intersectElement / unionElement;

	return jaccard;

}

	
double WMinHash::distance(WMinHash * wmh){
	return 1 - wJaccard(wmh);
}

WMinHash::~WMinHash()
{
	free(binsArr);
	free(countMinSketch);
	free(r);
	free(c);
	free(b);
	free(histoSketch_sketch);
	free(histoSketch_sketchWeight);

}

}