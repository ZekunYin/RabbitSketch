#ifndef __PYBIND_H__
#define __PYBIND_H__
#include <pybind11/pybind11.h>
namespace py = pybind11;

PYBIND11_MODULE(rabbitsketch, m){
	m.doc() = "rabbitsketch pybind";
	//py::class_<Sketch::Parameters>(m, "Parameters")
	//	.def(py::init<Sketch::Parameters &>())
	//	.def(py::init<>())
	//	.def_readwrite("kmerSize", &Sketch::Parameters::kmerSize)
	//	.def_readwrite("alphabetSize", &Sketch::Parameters::alphabetSize)
	//	.def_readwrite("preserveCase", &Sketch::Parameters::preserveCase)
	//	.def_readwrite("use64", &Sketch::Parameters::use64)
	//	.def_readwrite("sketchSize", &Sketch::Parameters::minHashesPerWindow)
	//	.def_readwrite("noncanonical", &Sketch::Parameters::noncanonical)
	//	//.def_readwrite("bumBins", &Sketch::Parameters::numBins)
	//	//.def_readwrite("minimizerWindowSize", &Sketch::Parameters::minimizerWindowSize)
	//	//.def_readwrite("histoSketch_sketchSize", &Sketch::Parameters::histoSketch_sketchSize)
	//	//.def_readwrite("histoSketch_dimension", &Sketch::Parameters::histoSketch_dimension)
	//	//.def_readwrite("paraDecayWeight", &Sketch::Parameters::paraDecayWeight)
	//	.def_readwrite("l", &Sketch::Parameters::l)
	//	.def_readwrite("m", &Sketch::Parameters::m)
	//	.def_readwrite("rc", &Sketch::Parameters::rc)
	//	;

	py::class_<Sketch::MinHash>(m, "MinHash")
		//.def(py::init<>())
		//.def(py::init<int>(), py::arg("k"))
		//.def(py::init<int, int>(), py::arg("k"), py::kwonly(), py::arg("size"))
		.def(py::init<int, int, uint32_t>(), py::arg("kmer") = 21, py::arg("size")=1000, py::arg("seed")=42)
		//.def(py::init<int, int, uint32_t>())
		.def("update", &Sketch::MinHash::update)
		.def("merge", &Sketch::MinHash::merge)
		.def("jaccard", &Sketch::MinHash::jaccard)
		.def("mdistance", &Sketch::MinHash::mdistance)
		.def("getTotalLength", &Sketch::MinHash::getTotalLength)
		.def("printMinHashes", &Sketch::MinHash::printMinHashes)
		.def("count", &Sketch::MinHash::count)

		//parameters
		//.def("setKmerSize", &Sketch::MinHash::setKmerSize)
		//.def("setAlphabetSize", &Sketch::MinHash::setAlphabetSize)
		//.def("setPreserveCase", &Sketch::MinHash::setPreserveCase)
		//.def("setUse64", &Sketch::MinHash::setUse64)
		//.def("setSeed", &Sketch::MinHash::setSeed)
		//.def("setSketchSize", &Sketch::MinHash::setSketchSize)
		//.def("setNoncanonical", &Sketch::MinHash::setNoncanonical)
		.def("getKmerSize", &Sketch::MinHash::getKmerSize)
		.def("getSeed", &Sketch::MinHash::getSeed)
		.def("getMaxSketchSize", &Sketch::MinHash::getMaxSketchSize)
		.def("getSketchSize", &Sketch::MinHash::getSketchSize)
		.def("isEmpty", &Sketch::MinHash::isEmpty)
		.def("isReverseComplement", &Sketch::MinHash::isReverseComplement)
		;

	py::class_<Sketch::WMinHash>(m, "WMinHash")
		//.def(py::init<>())
		.def(py::init<int, int, int, double>(), py::arg("kmer")=21, py::arg("size")=50, py::arg("windowSize")=9, py::arg("paraDWeight")=0.0)
		.def("update", &Sketch::WMinHash::update)
		.def("wJaccard", &Sketch::WMinHash::wJaccard)
		.def("distance", &Sketch::WMinHash::distance)
		.def("getWMinHash", &Sketch::WMinHash::getWMinHash)
		;

	py::class_<Sketch::OrderMinHash>(m, "OrderMinHash")
		.def(py::init<>())
		.def(py::init<char *>())
		.def("buildSketch", &Sketch::OrderMinHash::buildSketch)
		.def("similarity", &Sketch::OrderMinHash::similarity)
		.def("distance", &Sketch::OrderMinHash::distance)
		// parameters
		.def("setK", &Sketch::OrderMinHash::setK) 
		.def("setL", &Sketch::OrderMinHash::setL) 
		.def("setM", &Sketch::OrderMinHash::setM) 
		.def("setSeed", &Sketch::OrderMinHash::setSeed)
		.def("setReverseComplement", &Sketch::OrderMinHash::setReverseComplement)
		.def("getK", &Sketch::OrderMinHash::getK)
		.def("getL", &Sketch::OrderMinHash::getL)
		.def("getM", &Sketch::OrderMinHash::getM)
		.def("getSeed", &Sketch::OrderMinHash::getSeed)
		.def("isReverseComplement", &Sketch::OrderMinHash::isReverseComplement)
		;
}

#endif //__PYBIND_H__
