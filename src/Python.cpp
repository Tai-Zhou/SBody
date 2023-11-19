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

int InterpolatePosition(double result[], double position0[], const double position1[], double time_diff0, double time_diff1) {
#ifndef GSL_RANGE_CHECK_OFF
	py::print("last_position:", position0[0], position0[1], position0[2], position0[3], position0[4], position0[5], position0[6], position0[7]);
	py::print("this_position:", position1[0], position1[1], position1[2], position1[3], position1[4], position1[5], position1[6], position1[7]);
#endif
	const double coeff = 1. / (time_diff0 + time_diff1);
	copy(position1, position1 + 8, result);
	SphericalToCartesian(position0);
	SphericalToCartesian(result);
#ifndef GSL_RANGE_CHECK_OFF
	py::print("last_position:", position0[0], position0[1], position0[2], position0[3], position0[4], position0[5], position0[6], position0[7]);
	py::print("this_position:", result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7]);
#endif
	for (int i = 0; i < 8; ++i)
		result[i] = (time_diff1 * position0[i] + time_diff0 * result[i]) * coeff;
	return GSL_SUCCESS;
}

py::array_t<double>
CalculateFullStarOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e, double inclination, double ascending_node, double periapsis) {
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
	} else {
		main_metric = make_shared<Schwarzschild>();
		viewPtr = make_unique<View>(make_unique<Schwarzschild>(), R, 0.);
	}
	Star star_0(main_metric, T, LAGRANGIAN, GEODESIC, Unit::R_sun, false);
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
	auto result = py::array_t<double>(tStepNumber * 14);
	double *result_ptr = const_cast<double *>(result.data());
	for (int i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		status = star_0.IntegratorApply(&t, tStep, &h);
		if (status > 0)
			PrintlnError("main status = {}", status);
		star_0.Position(result_ptr);
		if (metric)
			viewPtr->TraceStar(result_ptr, i, result_ptr + 10, false);
		SphericalToCartesian(result_ptr);
		result_ptr[8] = t / Unit::s;
		const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
		result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
		result_ptr += 14;
	}
	return result.reshape({tStepNumber, 14UL});
}

py::array_t<double> CalculateStarOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e, double inclination, double ascending_node, double periapsis, bool ray_tracing, bool gr_time_delay, const py::array_t<double> &obs_time) {
#ifndef GSL_RANGE_CHECK_OFF
	py::print(mass, metric, fSP, R, tp, a, e, inclination, ascending_node, periapsis, ray_tracing, gr_time_delay);
#endif
	Unit::Initialize(mass);
	R *= Unit::pc;
	a *= R * Unit::mas;
	inclination *= M_PI / 180.;
	ascending_node *= M_PI / 180.;
	periapsis *= M_PI / 180.;
	double t = (tp - 2002.) * Unit::yr - M_PI * sqrt(gsl_pow_3(a)), tStep = 0., tRec = 0.001 * Unit::yr;
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
		viewPtr->TraceStar(position, 0, position + 8, false);
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
	double *last_position = new double[12], *this_position = new double[12];
	for (int i = 0;; ++i) {
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
				viewPtr->TraceStar(last_position, i, last_position + 8, false);
				viewPtr->TraceStar(this_position, i, this_position + 8, false);
#ifndef GSL_RANGE_CHECK_OFF
				py::print("gr_redshift: ", idx, last_position[10], this_position[10]);
#endif
				double this_gr_obs_time = (t + (t0 - this_position[11]) * Unit::s) / Unit::yr + 2002.;
				double last_gr_obs_time = (t - tRec + (t0 - last_position[11]) * Unit::s) / Unit::yr + 2002.;
#ifndef GSL_RANGE_CHECK_OFF
				py::print("gr_obs_time: ", idx, last_gr_obs_time, obs_time.at(idx), this_gr_obs_time);
#endif
				const double time_diff0 = obs_time.at(idx) - last_gr_obs_time, time_diff1 = this_gr_obs_time - obs_time.at(idx), coeff = 1. / (this_gr_obs_time - last_gr_obs_time);
#ifndef GSL_RANGE_CHECK_OFF
				py::print(last_gr_obs_time, obs_time.at(idx), this_gr_obs_time, time_diff1, time_diff0);
#endif
				InterpolatePosition(result_ptr, last_position, this_position, time_diff0, time_diff1);
				result_ptr[8] = (tStep - tRec * time_diff1 * coeff) / Unit::s;
				for (int j = 0; j < 4; ++j)
					result_ptr[10 + j] = (time_diff1 * last_position[j + 8] + time_diff0 * this_position[j + 8]) * coeff;
			} else {
				const double time_diff0 = obs_time.at(idx) - last_obs_time, time_diff1 = this_obs_time - obs_time.at(idx), coeff = 1. / (this_obs_time - last_obs_time);
				InterpolatePosition(result_ptr, last_position, this_position, time_diff0, time_diff1);
				if (ray_tracing || metric == 2) {
					CartesianToSpherical(result_ptr);
					if (ray_tracing)
						viewPtr->TraceStar(result_ptr, i, result_ptr + 10, false);
					if (metric == 2)
						result_ptr[1] -= 1.;
					SphericalToCartesian(result_ptr);
				}
				result_ptr[8] = (tStep - tRec * time_diff1 * coeff) / Unit::s;
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

double StarChi2(py::array_t<double> x, int metric, int gr_switch, int fix_switch, py::array_t<double> obs_time, py::array_t<double> obs_redshift, double redshift_sigma, py::array_t<double> obs_ra, double ra_sigma, py::array_t<double> obs_dec, double dec_sigma) {
	if (x.at(9) > 0.99)
		return GSL_NEGINF;
	if (metric == 0 && gr_switch > 0) {
		x.mutable_at(7) = 1.;
		gr_switch = 0;
	}
	if (fix_switch & 1)
		x.mutable_at(1) = 0.07;
	if (fix_switch & 2)
		x.mutable_at(2) = -0.9;
	if (fix_switch & 4)
		x.mutable_at(4) = 0.0341;
	if (fix_switch & 8)
		x.mutable_at(5) = 0.080;
	if (fix_switch & 16)
		x.mutable_at(6) = -1.6;
	auto obs_data = CalculateStarOrbit(x.at(0), metric, x.at(7), x.at(3), x.at(13), x.at(8), x.at(9), x.at(10), x.at(11), x.at(12), gr_switch > 0, (gr_switch & 4) > 0, obs_time);
	double redshift_prob = 0., ra_prob = 0., dec_prob = 0.;
	const int size = obs_redshift.size();
	if (gr_switch & 1)
		for (int i = 0; i < size; ++i) {
			redshift_prob += gsl_pow_2(obs_redshift.at(i) - (obs_data.at(i, 12) - 1.) * 299792.458 - x.at(6));
#ifndef GSL_RANGE_CHECK_OFF
			py::print(obs_redshift.at(i), (obs_data.at(i, 12) - 1.) * 299792.458 - x.at(6));
#endif
		}
	else
		for (int i = 0; i < size; ++i)
			redshift_prob += gsl_pow_2(obs_redshift.at(i) - obs_data.at(i, 9) - x.at(6));
	if (gr_switch & 2)
		for (int i = 0; i < size; ++i) {
			ra_prob += gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 10) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - (obs_time.at(i) - 2009.02) * x.at(5));
			dec_prob += gsl_pow_2(obs_dec.at(i) - obs_data.at(i, 11) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - (obs_time.at(i) - 2009.02) * x.at(4));
		}
	else
		for (int i = 0; i < size; ++i) {
			ra_prob += gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - (obs_time.at(i) - 2009.02) * x.at(5));
			dec_prob += gsl_pow_2(obs_dec.at(i) + obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - (obs_time.at(i) - 2009.02) * x.at(4));
#ifndef GSL_RANGE_CHECK_OFF
			py::print(obs_ra.at(i), obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - (obs_time.at(i) - 2009.02) * x.at(5));
			py::print(obs_dec.at(i), obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - (obs_time.at(i) - 2009.02) * x.at(4));
#endif
		}
#ifndef GSL_RANGE_CHECK_OFF
	py::print(">>>", redshift_prob / gsl_pow_2(redshift_sigma), ra_prob / gsl_pow_2(ra_sigma), dec_prob / gsl_pow_2(dec_sigma));
#endif
	double chi2 = redshift_prob / gsl_pow_2(redshift_sigma) + ra_prob / gsl_pow_2(ra_sigma) + dec_prob / gsl_pow_2(dec_sigma);
	if (isnan(chi2))
		for (int i = 0; i < 14; ++i)
			py::print(x.at(i));
	return chi2;
}

py::array_t<double> CalculateHSOrbit(double mass, int metric, int mode, double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination, bool ray_tracing, bool gr_time_delay, const py::array_t<double> &obs_time) {
	Unit::Initialize(mass);
	R *= Unit::pc;
	r *= R * Unit::mas;
	theta *= M_PI / 180.;
	phi *= M_PI / 180.;
	inclination *= M_PI / 180.;
	size_t tStepNumber = 10000;
	shared_ptr<Metric> main_metric;
	unique_ptr<View> viewPtr;
	if (metric == 0) {
		main_metric = make_shared<PN1>(fSP);
		ray_tracing = false;
		gr_time_delay = false;
	} else {
		main_metric = make_shared<Schwarzschild>();
		viewPtr = make_unique<View>(make_unique<Schwarzschild>(), R, inclination);
	}
	Star star_0(main_metric, T, LAGRANGIAN, GEODESIC, Unit::R_sun, false);
	double position[12], t = 0., tStep = 0., tRec = 1. * Unit::hr / tStepNumber, t0, last_obs_time = 0., this_obs_time;
	const double sin_inc = sin(inclination), cos_inc = cos(inclination);
	if (mode == 0)
		star_0.InitializeHelical(r, theta, phi, v_r, v_phi);
	else
		PrintlnWarning("Mode Error!");
	star_0.Position(position);
	double x0 = position[1] * abs(sin(position[2])) * cos(position[3]), z0 = position[1] * GSL_SIGN(position[2]) * cos(position[2]);
	if (gr_time_delay) {
		viewPtr->TraceStar(position, 0, position + 8, false);
		t0 = position[11];
	}
	double vz0 = x0 * sin_inc + z0 * cos_inc;
	double h = -1.;
	int status = 0;
	if (status = star_0.IntegratorApply(&t, 0., &h); status != 0)
		py::print("[!] IntegratorApply status =", status);
	h = 1.;
	if (status = star_0.IntegratorReset(); status != 0)
		py::print("[!] IntegratorReset status =", status);
	size_t idx = 0, size = obs_time.size();
	auto result = py::array_t<double>(size * 15);
	double *result_ptr = const_cast<double *>(result.data());
	double *last_position = new double[12], *this_position = new double[12], interpolation_position[8];
	for (int i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		if (status = star_0.IntegratorApply(&t, tStep, &h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		star_0.Position(this_position);
		this_obs_time = (t + vz0 - this_position[1] * (abs(sin(position[2])) * cos(position[3]) * sin_inc + GSL_SIGN(position[2]) * cos(position[2]) * cos_inc)) / Unit::yr + 2002.;
		if (this_obs_time > obs_time.at(idx)) {
			if (gr_time_delay) {
				viewPtr->TraceStar(last_position, i, last_position + 8, false);
				viewPtr->TraceStar(this_position, i, this_position + 8, false);
				double this_gr_obs_time = (t + (t0 - this_position[11]) * Unit::s) / Unit::yr + 2002.;
				double last_gr_obs_time = (t - tRec + (t0 - last_position[11]) * Unit::s) / Unit::yr + 2002.;
				const double time_diff0 = obs_time.at(idx) - last_gr_obs_time, time_diff1 = this_gr_obs_time - obs_time.at(idx), coeff = 1. / (this_gr_obs_time - last_gr_obs_time);
				InterpolatePosition(result_ptr, last_position, this_position, time_diff0, time_diff1);
				result_ptr[8] = (tStep - tRec * time_diff1 * coeff) / Unit::s;
				for (int j = 0; j < 4; ++j)
					result_ptr[10 + j] = (time_diff1 * last_position[j + 8] + time_diff0 * this_position[j + 8]) * coeff;
			} else {
				const double time_diff0 = obs_time.at(idx) - last_obs_time, time_diff1 = this_obs_time - obs_time.at(idx), coeff = 1. / (this_obs_time - last_obs_time);
				InterpolatePosition(result_ptr, last_position, this_position, time_diff0, time_diff1);
				if (ray_tracing || metric == 2) {
					CartesianToSpherical(result_ptr);
					if (ray_tracing)
						viewPtr->TraceStar(result_ptr, i, result_ptr + 10, false);
					if (metric == 2)
						result_ptr[1] -= 1.;
					SphericalToCartesian(result_ptr);
				}
				result_ptr[8] = (tStep - tRec * time_diff1 * coeff) / Unit::s;
			}
			const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
			result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
			result_ptr += 14;
			if (++idx >= size)
				break;
		}
		last_obs_time = this_obs_time;
		swap(last_position, this_position);
	}
	return result.reshape({size, 15UL});
}

double HSChi2(py::array_t<double> x, int metric, int mode, int gr_switch, int fix_switch, py::array_t<double> obs_time, py::array_t<double> obs_flux, double flux_sigma, py::array_t<double> obs_ra, double ra_sigma, py::array_t<double> obs_dec, double dec_sigma) {
	if (x.at(9) > 0.99)
		return GSL_NEGINF;
	if (metric == 0 && gr_switch > 0) {
		x.mutable_at(7) = 1.;
		gr_switch = 0;
	}
	auto obs_data = CalculateHSOrbit(x.at(0), metric, mode, x.at(7), x.at(3), x.at(8), x.at(9), x.at(10), x.at(11), x.at(12), x.at(13), gr_switch > 0, (gr_switch & 4) > 0, obs_time);
	double flux_prob = 0., ra_prob = 0., dec_prob = 0.;
	const int size = obs_flux.size();
	if (gr_switch & 1)
		for (int i = 0; i < size; ++i) {
			flux_prob += gsl_pow_2(obs_flux.at(i) - (obs_data.at(i, 12) - 1.) * 299792.458 - x.at(6));
		}
	else
		for (int i = 0; i < size; ++i)
			flux_prob += gsl_pow_2(obs_flux.at(i) - obs_data.at(i, 9) - x.at(6));
	if (gr_switch & 2)
		for (int i = 0; i < size; ++i) {
			ra_prob += gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 10) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - (obs_time.at(i) - 2009.02) * x.at(5));
			dec_prob += gsl_pow_2(obs_dec.at(i) - obs_data.at(i, 11) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - (obs_time.at(i) - 2009.02) * x.at(4));
		}
	else
		for (int i = 0; i < size; ++i) {
			ra_prob += gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - (obs_time.at(i) - 2009.02) * x.at(5));
			dec_prob += gsl_pow_2(obs_dec.at(i) + obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - (obs_time.at(i) - 2009.02) * x.at(4));
#ifndef GSL_RANGE_CHECK_OFF
			py::print(obs_ra.at(i), obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - (obs_time.at(i) - 2009.02) * x.at(5));
			py::print(obs_dec.at(i), obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - (obs_time.at(i) - 2009.02) * x.at(4));
#endif
		}
	double chi2 = flux_prob / gsl_pow_2(flux_sigma) + ra_prob / gsl_pow_2(ra_sigma) + dec_prob / gsl_pow_2(dec_sigma);
	if (isnan(chi2))
		for (int i = 0; i < 14; ++i)
			py::print(x.at(i));
	return chi2;
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

	m.def("CalculateFullStarOrbit", &CalculateFullStarOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("CalculateStarOrbit", &CalculateStarOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("StarChi2", &StarChi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

	// m.def("CalculateFullHSOrbit", &CalculateFull HSOrbit, R"pbdoc(
	//	Get trajectory, redshift, apparent position with orbital parameters
	// )pbdoc");

	m.def("CalculateHSOrbit", &CalculateHSOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("HSChi2", &HSChi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

	m.attr("s") = Unit::s;
	m.attr("day") = Unit::day;
	m.attr("yr") = Unit::yr;
	m.attr("__version__") = VERSION;
}
