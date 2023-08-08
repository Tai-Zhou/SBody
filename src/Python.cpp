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
#include <csignal>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include <fmt/core.h>
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

py::array_t<double> CalculateOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e, double inclination, double ascending_node, double periapsis, bool ray_tracing, py::array_t<double> obs_time) {
	Unit::Initialize(mass);
	R *= Unit::pc;
	a *= R * Unit::mas;
	inclination *= M_PI / 180.;
	ascending_node *= M_PI / 180.;
	periapsis *= M_PI / 180.;
	size_t tStepNumber = 10000;
	double t = (tp - 2002.) * Unit::yr - M_PI * sqrt(gsl_pow_3(a)), tStep = 0., tRec = 18. * Unit::yr / tStepNumber;
	shared_ptr<Metric> main_metric;
	unique_ptr<View> viewPtr;
	if (metric == 0)
		main_metric = make_shared<PN1>(fSP);
	else
		main_metric = make_shared<Schwarzschild>();
	Star star_0(main_metric, BASE, GEODESIC, Unit::R_sun, false);
	double position[8], z0, last_obs_time = 0., this_obs_time;
	star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, 0., 0., 0.);
	star_0.Position(position);
	z0 = position[3];
	star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
	if (metric != 0)
		viewPtr = make_unique<View>(make_unique<Schwarzschild>(), R, 0., "");
	double h = -1.;
	int status = 0;
	if (status = star_0.IntegratorApply(&t, 0., &h); status != 0)
		py::print("[!] IntegratorApply status =", status);
	h = 1.;
	if (status = star_0.IntegratorReset(); status != 0)
		py::print("[!] IntegratorReset status =", status);
	size_t idx = 0, size = obs_time.size();
	auto result = py::array_t<double>(size * 14);
	double *result_ptr = const_cast<double *>(result.data());
	double *last_position = new double[8], *this_position = new double[8];
	for (int i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		if (status = star_0.IntegratorApply(&t, tStep, &h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		star_0.Position(this_position);
		if (this_obs_time = (t + z0 - position[3]) / Unit::s; this_obs_time > obs_time.at(idx)) {
			auto interpolation_position = new double[8];
			const double coeff1 = this_obs_time - obs_time.at(idx), coeff2 = obs_time.at(idx) - last_obs_time, coeff3 = 1. / (this_obs_time - last_obs_time);
			for (int j = 0; j < 8; ++j)
				interpolation_position[j] = (coeff1 * last_position[j] + coeff2 * this_position[j]) * coeff3;
			if (ray_tracing)
				viewPtr->TraceStar(star_0, interpolation_position, i, result_ptr + 10);
			copy(interpolation_position, interpolation_position + 8, result_ptr);
			SphericalToCartesian(result_ptr);
			result_ptr[8] = (tStep - tRec * coeff1 * coeff3) / Unit::s;
			const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
			result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1) * 299792.458;
			result_ptr += 14;
			if (++idx >= size)
				break;
		}
		last_obs_time = this_obs_time;
		swap(last_position, this_position);
	}
	return result.reshape({size, 14UL});
}

double log_prob(double mass, int metric, double fSP, double R, double tp, double a, double e, double inclination, double ascending_node, double periapsis, bool gr_redshift, bool light_bending, bool gr_time_delay, py::array_t<double> obs_time, py::array_t<double> obs_redshift, double redshift_sigma, py::array_t<double> obs_ra, double ra_sigma, py::array_t<double> obs_dec, double dec_sigma) {
	auto obs_data = CalculateOrbit(mass, metric, fSP, R, tp, a, e, inclination, ascending_node, periapsis, gr_redshift || light_bending || gr_time_delay, obs_time);
	double redshift_prob = 0., ra_prob = 0., dec_prob = 0.;
	const int size = obs_redshift.size();
	if (gr_redshift)
		for (int i = 0; i < size; ++i)
			redshift_prob -= gsl_pow_2(obs_redshift.at(i) - obs_data.at(i, 12));
	else
		for (int i = 0; i < size; ++i)
			redshift_prob -= gsl_pow_2(obs_redshift.at(i) - obs_data.at(i, 9));
	if (light_bending)
		for (int i = 0; i < size; ++i) {
			ra_prob -= gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 10));
			dec_prob -= gsl_pow_2(obs_dec.at(i) - obs_data.at(i, 11));
		}
	else
		for (int i = 0; i < size; ++i) {
			ra_prob -= gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 2));
			dec_prob -= gsl_pow_2(obs_dec.at(i) + obs_data.at(i, 1));
		}
	return redshift_prob / gsl_pow_2(redshift_sigma) + ra_prob / gsl_pow_2(ra_sigma) + dec_prob / gsl_pow_2(dec_sigma);
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

	m.def("CalculateOrbit", &CalculateOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("log_prob", &log_prob, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

	m.attr("s") = Unit::s;
	m.attr("day") = Unit::day;
	m.attr("yr") = Unit::yr;
	m.attr("__version__") = VERSION;
}
