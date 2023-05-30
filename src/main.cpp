#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "ndarray.h"
#include "isosplit6.h"
#include "isocut6.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

int isosplit6_fn(py::array_t<int> labels_out, py::array_t<double> X, py::array_t<double> initial_labels)
{
    NDArray<double> Xa(X);
    NDArray<int32_t> initial_labels_a(labels_out);
    NDArray<int> La(labels_out);
    bigint N = Xa.shape[0];
    bigint M = Xa.shape[1];
    isosplit6_opts opts;
    isosplit6(La.ptr, M, N, Xa.ptr, initial_labels_a.ptr, opts);
    return 0;
}

int isocut6_fn(py::array_t<double> outputs, py::array_t<double> X)
{
    NDArray<double> Xa(X);
    NDArray<double> outputsa(outputs);
    bigint N = Xa.shape[0];
    isocut6_opts opts;

    double dipscore;
    double cutpoint;
    isocut6(&dipscore, &cutpoint, N, Xa.ptr, opts);
    outputsa.ptr[0] = dipscore;
    outputsa.ptr[1] = cutpoint;
    return 0;
}

namespace py = pybind11;

PYBIND11_MODULE(isosplit6_cpp, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: isosplit6_cpp

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("isosplit6_fn", &isosplit6_fn, "Isosplit6 clustering C++ implementation.");

    m.def("isocust6_fn", &isocut6_fn, "Isocut6 C++ implementation.");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
