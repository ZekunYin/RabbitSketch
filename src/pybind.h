#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(sketch, m){
	m.doc() = "sketch pybind";
	py::class_<Sketch::Parameters>(m, "Parameters")
		.def(py::init<Sketch::Parameters &>())
		.def(py::init<>())
		.def_readwrite("kmerSize", &Sketch::Parameters::kmerSize)
		.def_readwrite("alphabetSize", &Sketch::Parameters::alphabetSize)
		.def_readwrite("preserveCase", &Sketch::Parameters::preserveCase)
		.def_readwrite("use64", &Sketch::Parameters::use64)
		.def_readwrite("sketchSize", &Sketch::Parameters::minHashesPerWindow)
		.def_readwrite("noncanonical", &Sketch::Parameters::noncanonical)
		;
	py::class_<Sketch::MinHash>(m, "MinHash")
		.def(py::init<Sketch::Parameters>())
		.def("update", &Sketch::MinHash::update)
		.def("merge", &Sketch::MinHash::merge)
		.def("jaccard", &Sketch::MinHash::jaccard)
		//.def("dist", &Sketch::MinHash::dist)
		.def("getTotalLength", &Sketch::MinHash::getTotalLength)
		.def("printMinHashes", &Sketch::MinHash::printMinHashes)
		;
}
