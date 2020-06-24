#include "Sketch.h"
#include <random>
#include "histoSketch.h"

namespace Sketch{

//start the WMinHash@xxm
//WMinHash::WMinHash(Parameters parametersNew)//:parameters(parametersNew)
WMinHash::WMinHash()
{	
	//parameters = parametersNew;

	numBins = pow(kmerSize, alphabetSize); //need to be confirmed
	//histoDimension = pow(kmerSize, alphabetSize); //need to be confirmed

	//binsArr = (double *) malloc (parameters.get_numBins() * sizeof(double));
	binsArr = (double *)malloc(numBins * sizeof(double));
	for(int i = 0; i < numBins; i++){
		binsArr[i] = 0.0;
	}

	int g = ceil(2 / EPSILON);
	int d = ceil(log(1 - DELTA) / log(0.5));
	countMinSketch = (double *)malloc(d * g *sizeof(double));
	for(int i = 0; i < d * g; i++){
		countMinSketch[i] = 0.0;
	}
	
	//the r, c, b and getCWS need to be outClass
	r = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	c = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	b = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	getCWS(r, c, b, histoSketchSize, histoDimension);
	

	histoSketches = (uint32_t *) malloc (histoSketchSize * sizeof(uint32_t));
	histoWeight = (double *) malloc (histoSketchSize * sizeof(double));
	
	//add the applyConceptDrift and decayWeight.
	if(paraDecayWeight < 0.0 || paraDecayWeight > 1.0){
		cerr << "the paraDecayWeight must between 0.0 and 1.0 " << endl;
		exit(1);
	}
	else{
		applyConceptDrift = true;
	}
	if(paraDecayWeight == 1.0){
		applyConceptDrift = false;
	}
	decayWeight = 1.0;
	if(applyConceptDrift){
		decayWeight = exp(-paraDecayWeight);
	}

	needToCompute = true;

}

void WMinHash::update(char * seq)
{
	//int k = parameters.kmerSize;
	//int w = parameters.get_minimizerWindowSize();
	//int numBins = parameters.get_numBins();
	//int histoSketch_sketchSize = parameters.get_histoSketch_sketchSize();
	//int histoSketch_dimension = parameters.get_histoSketch_dimension();
	//cerr << "in the update the histoSketch_sketchSize is: " << parameters.get_histoSketch_sketchSize() << endl;

	//findMinimizers(k, w, seq, sketches);
	//kmerSpectrumAddHash(sketches, binsArr, parameters.get_numBins());
	findMinimizers(kmerSize, minimizerWindowSize, seq, sketches);
	kmerSpectrumAddHash(sketches, binsArr, numBins);
//	for(int i = 0; i < numBins; i++){
//		if(binsArr[i] != 0){
//			printf("%d\t%lf\n", i, binsArr[i]);
//		}
//	}
	sketches.clear();
	needToCompute = true;

}

void WMinHash::computeHistoSketch()
{
	kmerSpectrumDump(binsArr, numBins, kmerSpectrums);
	//cerr << "finish the kmerSpectrumDump " << endl;
	//int hsSketchSize = parameters.get_histoSketch_sketchSize();
	//int hsDimension = parameters.get_histoSketch_dimension();
	//cerr << "the hsSketchSize is " << hsSketchSize << endl;
	//cerr << "the hsDimension is " << hsDimension << endl;

	for(int i = 0; i < kmerSpectrums.size(); i++){
	//cerr << "finish the " << i << " iterator to addElement" << endl;
		//histoSketchAddElement((uint64_t)kmerSpectrums[i].BinID, kmerSpectrums[i].Frequency, countMinSketch, parameters.get_histoSketch_sketchSize(), applyConceptDrift, decayWeight, r, c, b, parameters.get_histoSketch_sketchSize(), parameters.get_histoSketch_dimension(), histoSketch_sketch, histoSketch_sketchWeight);
		histoSketchAddElement((uint64_t)kmerSpectrums[i].BinID, kmerSpectrums[i].Frequency, countMinSketch, applyConceptDrift, decayWeight, r, c, b, histoSketchSize, histoDimension, histoSketches, histoWeight);
	}
	//cerr << "finish the histoSketchAddElement " << endl;

}

void WMinHash::getWMinHash(){
	if(needToCompute){
		computeHistoSketch();
		needToCompute = false;
	}
	cout << "the sketch is: " << endl;
	for(int i = 0; i < histoSketchSize; i++){
		printf("%d\n", histoSketches[i]);
	}

	cout << "the sketchweight is: " << endl;
	for(int j = 0; j < histoSketchSize; j++){
		printf("%lf\n", histoWeight[j]);
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
	cerr << "the parameters.get_histoSketch_sketchSize() is: " << histoSketchSize << endl;
	for(int i = 0 ; i < histoSketchSize; i++){

		double curWeightA = abs(histoWeight[i]);
		double curWeightB = abs(wmh->histoWeight[i]);

		//get the intersection and union values
		if(histoSketches[i] == wmh->histoSketches[i]){
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

void WMinHash::setHistoSketchSize(int histoSketchSizeNew){
	histoSketchSize = histoSketchSizeNew;

	r = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	c = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	b = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	getCWS(r, c, b, histoSketchSize, histoDimension);
	

	histoSketches = (uint32_t *) malloc (histoSketchSize * sizeof(uint32_t));
	histoWeight = (double *) malloc (histoSketchSize * sizeof(double));
}

void WMinHash::setHistoDimension(int histoDimensionNew){
	histoDimension = histoDimensionNew;
	
	r = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	c = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	b = (double *)malloc(histoSketchSize * histoDimension * sizeof(double));
	getCWS(r, c, b, histoSketchSize, histoDimension);
	

	histoSketches = (uint32_t *) malloc (histoSketchSize * sizeof(uint32_t));
	histoWeight = (double *) malloc (histoSketchSize * sizeof(double));
}


WMinHash::~WMinHash()
{
	free(binsArr);
	free(countMinSketch);
	free(r);
	free(c);
	free(b);
	free(histoSketches);
	free(histoWeight);

}

}
