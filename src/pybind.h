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
		.def_readwrite("bumBins", &Sketch::Parameters::numBins)
		.def_readwrite("minimizerWindowSize", &Sketch::Parameters::minimizerWindowSize)
		.def_readwrite("histoSketch_sketchSize", &Sketch::Parameters::histoSketch_sketchSize)
		.def_readwrite("histoSketch_dimension", &Sketch::Parameters::histoSketch_dimension)
		.def_readwrite("paraDecayWeight", &Sketch::Parameters::paraDecayWeight)
		.def_readwrite("l", &Sketch::Parameters::l)
		.def_readwrite("m", &Sketch::Parameters::m)
		.def_readwrite("rc", &Sketch::Parameters::rc)
		;

	py::class_<Sketch::MinHash>(m, "MinHash")
		.def(py::init<Sketch::Parameters>())
		.def("update", &Sketch::MinHash::update)
		.def("merge", &Sketch::MinHash::merge)
		.def("jaccard", &Sketch::MinHash::jaccard)
		.def("dist", &Sketch::MinHash::dist)
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
		.def(py::init<Sketch::Parameters, char *>())
		.def("sketch", &Sketch::OMinHash::sketch)
		.def("similarity", &Sketch::OMinHash::similarity)
		.def("distance", &Sketch::OMinHash::distance)
		;
}
