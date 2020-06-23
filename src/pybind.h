#ifndef __PYBIND_H__
#define __PYBIND_H__
#include <pybind11/pybind11.h>
namespace py = pybind11;

PYBIND11_MODULE(rabbitsketch, m){
	m.doc() = "rabbitsketch pybind";
	py::class_<Sketch::Parameters>(m, "Parameters")
		.def(py::init<Sketch::Parameters &>())
		.def(py::init<>())
		.def_readwrite("kmerSize", &Sketch::Parameters::kmerSize)
		.def_readwrite("alphabetSize", &Sketch::Parameters::alphabetSize)
		.def_readwrite("preserveCase", &Sketch::Parameters::preserveCase)
		.def_readwrite("use64", &Sketch::Parameters::use64)
		.def_readwrite("sketchSize", &Sketch::Parameters::minHashesPerWindow)
		.def_readwrite("noncanonical", &Sketch::Parameters::noncanonical)
		//.def_readwrite("bumBins", &Sketch::Parameters::numBins)
		//.def_readwrite("minimizerWindowSize", &Sketch::Parameters::minimizerWindowSize)
		//.def_readwrite("histoSketch_sketchSize", &Sketch::Parameters::histoSketch_sketchSize)
		//.def_readwrite("histoSketch_dimension", &Sketch::Parameters::histoSketch_dimension)
		//.def_readwrite("paraDecayWeight", &Sketch::Parameters::paraDecayWeight)
		.def_readwrite("l", &Sketch::Parameters::l)
		.def_readwrite("m", &Sketch::Parameters::m)
		.def_readwrite("rc", &Sketch::Parameters::rc)
		;

	py::class_<Sketch::MinHash>(m, "MinHash")
		.def(py::init<Sketch::Parameters>())
		.def("update", &Sketch::MinHash::update)
		.def("merge", &Sketch::MinHash::merge)
		.def("jaccard", &Sketch::MinHash::jaccard)
		.def("distance", &Sketch::MinHash::distance)
		.def("getTotalLength", &Sketch::MinHash::getTotalLength)
		.def("printMinHashes", &Sketch::MinHash::printMinHashes)
		;

	py::class_<Sketch::WMinHash>(m, "WMinHash")
		.def(py::init<Sketch::Parameters>())
		.def("update", &Sketch::WMinHash::update)
		.def("wJaccard", &Sketch::WMinHash::wJaccard)
		.def("distance", &Sketch::WMinHash::distance)
		.def("getWMinHash", &Sketch::WMinHash::getWMinHash)
		;

	py::class_<Sketch::OMinHash>(m, "OMinHash")
		.def(py::init<>())
		.def(py::init<char *>())
		.def("buildSketch", &Sketch::OMinHash::buildSketch)
		.def("similarity", &Sketch::OMinHash::similarity)
		.def("distance", &Sketch::OMinHash::distance)
		// parameters
		.def("setK", &Sketch::OMinHash::setK) 
		.def("setL", &Sketch::OMinHash::setL) 
		.def("setM", &Sketch::OMinHash::setM) 
		.def("setSeed", &Sketch::OMinHash::setSeed)
		.def("setReverseComplement", &Sketch::OMinHash::setReverseComplement)
		.def("getK", &Sketch::OMinHash::getK)
		.def("getL", &Sketch::OMinHash::getL)
		.def("getM", &Sketch::OMinHash::getM)
		.def("getSeed", &Sketch::OMinHash::getSeed)
		.def("isReverseComplement", &Sketch::OMinHash::isReverseComplement)
		;
}

#endif //__PYBIND_H__