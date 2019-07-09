#include <iostream>
#include "Sketch.h"

using namespace std;

int main(){
	Sketch::Parameters parameter;
	
    parameter.kmerSize = 21;
	parameter.minHashesPerWindow = 1000;
    parameter.seed = 42;
	parameter.use64 = true;
	char seq[] = "AGCTAAACGGTACGATCTGTGCAT";
	Sketch::SketchInput * input = new Sketch::SketchInput(seq , 24, parameter);
	sketchOutput(input);
}
