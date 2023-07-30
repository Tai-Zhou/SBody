/**
 * @file Python.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <array>
#include <chrono>
#include <cmath>
#include <csignal>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include <gsl/gsl_math.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Unit.h"
#include "Utility.h"
#include "View.h"

namespace py = pybind11;
using namespace std;
using namespace SBody;

py::array_t<double> MCMC(double mass, int metric, double PN, double R, double tp, double a = 5., double e = 0.884649, double inclination = 134.567, double ascending_node = 228.171, double periapsis = 66.263) {
	Unit::Initialize(mass);
	a *= Unit::mas * R * Unit::pc;
	inclination *= M_PI / 180.;
	ascending_node *= M_PI / 180.;
	periapsis *= M_PI / 180.;
	double true_anomaly = M_PI;
	size_t tStepNumber = 10000;
	double t = (tp - 2002.) * Unit::yr - M_PI * sqrt(gsl_pow_3(a)), tStep = 0., tRec = 18. * Unit::yr / tStepNumber;
	shared_ptr<Metric> main_metric = make_shared<Newton>(1, T);
	Star star_0(main_metric, Unit::R_sun, 0);
	star_0.InitializeKeplerian(a * Unit::mpc, e, inclination, periapsis, ascending_node, true_anomaly, 0., 0.);
	auto result = py::array_t<double>(tStepNumber * 9);
	double *result_ptr = static_cast<double *>(result.request().ptr);
	double h = -1.;
	int status = 0;
	if (status = star_0.IntegratorApply(&t, 0., &h); status != 0)
		py::print("[!] IntegratorApply status =", status);
	h = 1.;
	if (status = star_0.IntegratorReset(); status != 0)
		py::print("[!] IntegratorReset status =", status);
	for (int i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		if (status = star_0.IntegratorApply(&t, tStep, &h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		star_0.Position(result_ptr);
		if (metric != 0)
			SphericalToCartesian(result_ptr);
		const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
		result_ptr[8] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1) * 299792.458;
		result_ptr += 9;
	}
	return result.reshape({tStepNumber, 9UL});
}

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

	m.def("MCMC", &MCMC, R"pbdoc(
		MCMC function
	)pbdoc");

	m.attr("s") = Unit::s;
	m.attr("day") = Unit::day;
	m.attr("yr") = Unit::yr;
	m.attr("__version__") = VERSION;
}
