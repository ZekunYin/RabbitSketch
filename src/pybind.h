#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(sketch, m){
	m.doc() = "sketch pybind";
	py::class_<Sketch::MinHash>(m, "MinHash")
		.def(py::init<Sketch::Parameters>())
		.def("jaccard", &Sketch::MinHash::jaccard);
}
