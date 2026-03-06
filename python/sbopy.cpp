/**
 * @file sbopy.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <pybind11/pybind11.h>

#include "infra/status.hpp"
#include "math/special_functions.hpp"

namespace py = pybind11;
using namespace sbody;

#define PRECISION double

void bindSStarModule(py::module_ &m);
void bindHotspotModule(py::module_ &m);

PYBIND11_MODULE(SBoPy, m) {
	m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: cmake_example
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

	m.attr("__version__") = VERSION;

	py::enum_<StatusCode>(m, "StatusCode")
		.value("SUCCESS", StatusCode::SUCCESS)
		.value("CONTINUE", StatusCode::CONTINUE)
		.value("FAILURE", StatusCode::FAILURE)
		.value("INVALID_ARGUMENT", StatusCode::INVALID_ARGUMENT)
		.value("DOMAIN_ERROR", StatusCode::DOMAIN_ERROR)
		.value("DIVISION_BY_ZERO", StatusCode::DIVISION_BY_ZERO)
		.value("MAX_ITERATIONS_EXCEEDED", StatusCode::MAX_ITERATIONS_EXCEEDED)
		.value("SINGULAR_MATRIX", StatusCode::SINGULAR_MATRIX)
		.value("NUMERIC_ERROR", StatusCode::NUMERIC_ERROR)
		.value("TRUST_REGION_TOO_SMALL", StatusCode::TRUST_REGION_TOO_SMALL)
		.value("SCALING_REQUIRED", StatusCode::SCALING_REQUIRED)
		.value("STEP_SIZE_TOO_SMALL", StatusCode::STEP_SIZE_TOO_SMALL)
		.value("NO_INTERSECTION_WITH_TARGET_SURFACE", StatusCode::NO_INTERSECTION_WITH_TARGET_SURFACE)
		.value("FALL_INTO_BLACK_HOLE", StatusCode::FALL_INTO_BLACK_HOLE)
		.value("INPUT_OUTPUT_ERROR", StatusCode::INPUT_OUTPUT_ERROR)
		.value("NOT_IMPLEMENTED", StatusCode::NOT_IMPLEMENTED)
		.value("COMPUTATION_TIMEOUT", StatusCode::COMPUTATION_TIMEOUT);

	bindSStarModule(m);
	bindHotspotModule(m);

	m.def("EllipticIntegral", &EllipticIntegral<double>, R"pbdoc(
		Elliptic Integral
	)pbdoc");

	m.def("EllipticIntegral2Complex", &EllipticIntegral2Complex<double>, R"pbdoc(
		Elliptic Integral with 2 complex roots
	)pbdoc");

	m.def("EllipticIntegral4Complex", &EllipticIntegral4Complex<double>, R"pbdoc(
		Elliptic Integral with 4 complex roots
	)pbdoc");
}
