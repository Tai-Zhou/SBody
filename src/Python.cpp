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

py::array_t<double> CalculateOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e, double inclination, double ascending_node, double periapsis, bool ray_tracing, bool gr_time_delay, const py::array_t<double> &obs_time) {
#ifndef GSL_RANGE_CHECK_OFF
	py::print(mass, metric, fSP, R, tp, a, e, inclination, ascending_node, periapsis, ray_tracing, gr_time_delay);
#endif
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
	if (metric == 0) {
		main_metric = make_shared<PN1>(fSP);
		ray_tracing = false;
		gr_time_delay = false;
	} else {
		main_metric = make_shared<Schwarzschild>();
		viewPtr = make_unique<View>(make_unique<Schwarzschild>(), R, 0.);
	}
	Star star_0(main_metric, T, LAGRANGIAN, GEODESIC, Unit::R_sun, false);
	double position[12], z0, t0, last_obs_time = 2002., this_obs_time;
	if (metric == 2)
		star_0.InitializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, 0., 0., 0.);
	else
		star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, 0., 0., 0.);
	star_0.Position(position);
#ifndef GSL_RANGE_CHECK_OFF
	py::print(position[0], position[1], position[2], position[3], position[4], position[5], position[6], position[7]);
#endif
	z0 = position[1] * GSL_SIGN(position[2]) * cos(position[2]);
	if (gr_time_delay) {
		viewPtr->TraceStar(star_0, nullptr, 0, position + 8);
		t0 = position[11];
	}
	if (metric == 2)
		star_0.InitializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
	else
		star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
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
	double *last_position = new double[12], *this_position = new double[12], interpolation_position[8];
	for (int i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		if (status = star_0.IntegratorApply(&t, tStep, &h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		star_0.Position(this_position);
		if (metric == 2)
			this_obs_time = (t + z0 - (this_position[1] - 1.) * GSL_SIGN(this_position[2]) * cos(this_position[2])) / Unit::yr + 2002.;
		else
			this_obs_time = (t + z0 - this_position[1] * GSL_SIGN(this_position[2]) * cos(this_position[2])) / Unit::yr + 2002.;
		if (this_obs_time > obs_time.at(idx)) {
#ifndef GSL_RANGE_CHECK_OFF
			py::print("obs_time: ", idx, last_obs_time, obs_time.at(idx), this_obs_time);
#endif
			if (gr_time_delay) {
				viewPtr->TraceStar(star_0, last_position, i, last_position + 8);
				viewPtr->TraceStar(star_0, this_position, i, this_position + 8);
#ifndef GSL_RANGE_CHECK_OFF
				py::print("gr_redshift: ", idx, last_position[10], this_position[10]);
#endif
				double this_gr_obs_time = (t + (t0 - this_position[11]) * Unit::s) / Unit::yr + 2002.;
				double last_gr_obs_time = (t - tRec + (t0 - last_position[11]) * Unit::s) / Unit::yr + 2002.;
#ifndef GSL_RANGE_CHECK_OFF
				py::print("gr_obs_time: ", idx, last_gr_obs_time, obs_time.at(idx), this_gr_obs_time);
#endif
				const double coeff1 = this_gr_obs_time - obs_time.at(idx), coeff2 = obs_time.at(idx) - last_gr_obs_time, coeff3 = 1. / (this_gr_obs_time - last_gr_obs_time);
#ifndef GSL_RANGE_CHECK_OFF
				py::print(last_gr_obs_time, obs_time.at(idx), this_gr_obs_time, coeff1, coeff2);
#endif
				for (int j = 0; j < 8; ++j)
					result_ptr[j] = (coeff1 * last_position[j] + coeff2 * this_position[j]) * coeff3;
				SphericalToCartesian(result_ptr);
				result_ptr[8] = (tStep - tRec * coeff1 * coeff3) / Unit::s;
				for (int j = 0; j < 4; ++j)
					result_ptr[10 + j] = (coeff1 * last_position[j + 8] + coeff2 * this_position[j + 8]) * coeff3;
			} else {
				const double coeff1 = this_obs_time - obs_time.at(idx), coeff2 = obs_time.at(idx) - last_obs_time, coeff3 = 1. / (this_obs_time - last_obs_time);
				for (int j = 0; j < 8; ++j)
					interpolation_position[j] = (coeff1 * last_position[j] + coeff2 * this_position[j]) * coeff3;
#ifndef GSL_RANGE_CHECK_OFF
				py::print("last_position:", last_position[0], last_position[1], last_position[2], last_position[3], last_position[4], last_position[5], last_position[6], last_position[7]);
				py::print("this_position:", this_position[0], this_position[1], this_position[2], this_position[3], this_position[4], this_position[5], this_position[6], this_position[7]);
				SphericalToCartesian(last_position);
				SphericalToCartesian(this_position);
				double interpolation_position2[12], interpolation_position3[8];
				for (int j = 0; j < 8; ++j)
					interpolation_position3[j] = interpolation_position2[j] = (coeff1 * last_position[j] + coeff2 * this_position[j]) * coeff3;
				CartesianToSpherical(this_position);
				CartesianToSpherical(interpolation_position3);
				viewPtr->TraceStar(star_0, interpolation_position3, i, interpolation_position2 + 8);
				py::print("interpolation_position2:", interpolation_position2[0], interpolation_position2[1], interpolation_position2[2], interpolation_position2[3], interpolation_position2[4], interpolation_position2[5], interpolation_position2[6], interpolation_position2[7], interpolation_position2[8], interpolation_position2[9], interpolation_position2[10], interpolation_position2[11]);
#endif
				if (ray_tracing)
					viewPtr->TraceStar(star_0, interpolation_position, i, result_ptr + 10);
				copy(interpolation_position, interpolation_position + 8, result_ptr);
				if (metric == 2)
					result_ptr[1] -= 1.;
				SphericalToCartesian(result_ptr);
				result_ptr[8] = (tStep - tRec * coeff1 * coeff3) / Unit::s;
			}
			const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
			result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
#ifndef GSL_RANGE_CHECK_OFF
			py::print(result_ptr[0], result_ptr[1], result_ptr[2], result_ptr[3], result_ptr[4], result_ptr[5], result_ptr[6], result_ptr[7], result_ptr[8], result_ptr[9], result_ptr[10], result_ptr[11], result_ptr[12], result_ptr[13]);
#endif
			result_ptr += 14;
			if (++idx >= size)
				break;
		}
		last_obs_time = this_obs_time;
		swap(last_position, this_position);
	}
	return result.reshape({size, 14UL});
}

double Chi2(py::array_t<double> x, int metric, int gr_switch, py::array_t<double> obs_time, py::array_t<double> obs_redshift, double redshift_sigma, py::array_t<double> obs_ra, double ra_sigma, py::array_t<double> obs_dec, double dec_sigma) {
	if (x.at(9) > 0.99)
		return GSL_NEGINF;
	if (metric == 0)
		gr_switch = 0;
	auto obs_data = CalculateOrbit(x.at(0), metric, x.at(7), x.at(3), x.at(13), x.at(8), x.at(9), x.at(10), x.at(11), x.at(12), gr_switch > 0, (gr_switch & 4) > 0, obs_time);
	double redshift_prob = 0., ra_prob = 0., dec_prob = 0.;
	const int size = obs_redshift.size();
	if (gr_switch & 1)
		for (int i = 0; i < size; ++i) {
			redshift_prob -= gsl_pow_2(obs_redshift.at(i) - (obs_data.at(i, 12) - 1.) * 299792.458 - x.at(6));
#ifndef GSL_RANGE_CHECK_OFF
			py::print(obs_redshift.at(i), (obs_data.at(i, 12) - 1.) * 299792.458 - x.at(6));
#endif
		}
	else
		for (int i = 0; i < size; ++i)
			redshift_prob -= gsl_pow_2(obs_redshift.at(i) - obs_data.at(i, 9) - x.at(6));
	if (gr_switch & 2)
		for (int i = 0; i < size; ++i) {
			ra_prob -= gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 10) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - obs_time.at(i) * x.at(5));
			dec_prob -= gsl_pow_2(obs_dec.at(i) - obs_data.at(i, 11) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - obs_time.at(i) * x.at(4));
		}
	else
		for (int i = 0; i < size; ++i) {
			ra_prob -= gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - obs_time.at(i) * x.at(5));
			dec_prob -= gsl_pow_2(obs_dec.at(i) + obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - obs_time.at(i) * x.at(4));
#ifndef GSL_RANGE_CHECK_OFF
			py::print(obs_ra.at(i), obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - obs_time.at(i) * x.at(5));
			py::print(obs_dec.at(i), obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - obs_time.at(i) * x.at(4));
#endif
		}
#ifndef GSL_RANGE_CHECK_OFF
	py::print(">>>", redshift_prob / gsl_pow_2(redshift_sigma), ra_prob / gsl_pow_2(ra_sigma), dec_prob / gsl_pow_2(dec_sigma));
#endif
	double ans = redshift_prob / gsl_pow_2(redshift_sigma) + ra_prob / gsl_pow_2(ra_sigma) + dec_prob / gsl_pow_2(dec_sigma);
	if (isnan(ans))
		for (int i = 0; i < 14; ++i)
			py::print(x.at(i));
	return ans;
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

	m.def("Chi2", &Chi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

	m.attr("s") = Unit::s;
	m.attr("day") = Unit::day;
	m.attr("yr") = Unit::yr;
	m.attr("__version__") = VERSION;
}
