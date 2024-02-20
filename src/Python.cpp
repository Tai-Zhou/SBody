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

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Unit.h"
#include "Utility.h"
#include "View.h"

namespace py = pybind11;
using namespace std;
using namespace SBody;

double CalculatePericenterTime(double mass, int metric, double fSP, double R, double ta, double a, double e, double inclination, double ascending_node, double periapsis) {
	Unit::Initialize(mass);
	R *= Unit::pc;
	a *= R * Unit::mas;
	inclination *= M_PI / 180.;
	ascending_node *= M_PI / 180.;
	periapsis *= M_PI / 180.;
	double t = (ta - 2002.) * Unit::yr, tStep = 0., tRec = 0.000001 * Unit::yr;
	shared_ptr<Metric> main_metric;
	if (metric == 0)
		main_metric = make_shared<PN1>(fSP);
	else
		main_metric = make_shared<Schwarzschild>();
	Particle star_0(main_metric, T, LAGRANGIAN, false);
	if (metric == 2)
		star_0.InitializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
	else
		star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
	double h = 1., position[8], last_radius = GSL_POSINF;
	star_0.Position(position);
	int status = 0;
	if (status = star_0.IntegratorApply(&t, 32.3 * Unit::yr, &h); status != 0)
		py::print("[!] IntegratorApply status =", status);
	for (int i = 0;; ++i) {
		status = star_0.IntegratorApply(&t, tStep, &h);
		if (status > 0)
			PrintlnError("main status = {}", status);
		star_0.Position(position);
		if (position[1] > last_radius)
			return tStep / Unit::yr + 2002.;
		last_radius = position[1];
		tStep += tRec;
	}
}

py::array_t<double> CalculateFullStarOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e, double inclination, double ascending_node, double periapsis, double t1) {
	Unit::Initialize(mass);
	R *= Unit::pc;
	a *= R * Unit::mas;
	inclination *= M_PI / 180.;
	ascending_node *= M_PI / 180.;
	periapsis *= M_PI / 180.;
	size_t tStepNumber = 1000;
	double t = (tp - 2002.) * Unit::yr - M_PI * sqrt(gsl_pow_3(a)), tStep = 0., tRec = (t1 - 2002.) * Unit::yr / tStepNumber;
	shared_ptr<Metric> main_metric;
	unique_ptr<View> view_ptr;
	if (metric == 0) {
		main_metric = make_shared<PN1>(fSP);
	} else {
		main_metric = make_shared<Schwarzschild>();
		view_ptr = make_unique<View>(make_unique<Schwarzschild>(), R, 0., 0.);
	}
	time_system star_time = T;
	Particle star_0(main_metric, star_time, LAGRANGIAN, false);
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
	for (size_t i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		status = star_0.IntegratorApply(&t, tStep, &h);
		if (status > 0)
			PrintlnError("main status = {}", status);
		star_0.Position(result_ptr);
		if (metric)
			view_ptr->Trace(result_ptr, star_time, result_ptr + 10, false);
		SphericalToCartesian(result_ptr);
		result_ptr[8] = t / Unit::s;
		const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
		result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
		result_ptr += 14;
	}
	return result.reshape({tStepNumber, 14UL});
}

py::array_t<double> CalculateStarOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e, double inclination, double ascending_node, double periapsis, bool ray_tracing, bool gr_time_delay, const py::array_t<double> &obs_time) {
	Unit::Initialize(mass);
	R *= Unit::pc;
	a *= R * Unit::mas;
	inclination *= M_PI / 180.;
	ascending_node *= M_PI / 180.;
	periapsis *= M_PI / 180.;
	double t = (tp - 2002.) * Unit::yr - M_PI * sqrt(gsl_pow_3(a)), tStep = 0., tRec = 0.001 * Unit::yr;
	shared_ptr<Metric> main_metric;
	unique_ptr<View> view_ptr;
	if (metric == 0) {
		main_metric = make_shared<PN1>(fSP);
		ray_tracing = false;
		gr_time_delay = false;
	} else {
		main_metric = make_shared<Schwarzschild>();
		view_ptr = make_unique<View>(make_unique<Schwarzschild>(), R, 0., 0.);
	}
	time_system star_time = T;
	Particle star_0(main_metric, star_time, LAGRANGIAN, false);
	double z0, t0, last_obs_time = 2002., this_obs_time, *last_position = new double[13], *this_position = new double[13];
	if (metric == 2)
		star_0.InitializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, 0., 0., 0.);
	else
		star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, 0., 0., 0.);
	star_0.Position(last_position);
	z0 = last_position[1] * GSL_SIGN(last_position[2]) * cos(last_position[2]);
	if (gr_time_delay) {
		view_ptr->Trace(last_position, star_time, last_position + 8, false);
		t0 = last_position[11];
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
			if (gr_time_delay) {
				view_ptr->Trace(last_position, star_time, last_position + 8, false);
				view_ptr->Trace(this_position, star_time, this_position + 8, false);
				double this_gr_obs_time = (t + (t0 - this_position[11]) * Unit::s) / Unit::yr + 2002.;
				double last_gr_obs_time = (t - tRec + (t0 - last_position[11]) * Unit::s) / Unit::yr + 2002.;
				InterpolateSphericalPositionToCartesian(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, last_position, this_position, result_ptr);
				result_ptr[8] = LinearInterpolation(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, tStep - tRec, tStep) / Unit::s;
				LinearInterpolation(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, last_position + 8, this_position + 8, result_ptr + 10, 4);
			} else {
				InterpolateSphericalPositionToCartesian(obs_time.at(idx), last_obs_time, this_obs_time, last_position, this_position, result_ptr);
				if (ray_tracing || metric == 2) {
					CartesianToSpherical(result_ptr);
					if (ray_tracing)
						view_ptr->Trace(result_ptr, star_time, result_ptr + 10, false);
					if (metric == 2)
						result_ptr[1] -= 1.;
					SphericalToCartesian(result_ptr);
				}
				result_ptr[8] = LinearInterpolation(obs_time.at(idx), last_obs_time, this_obs_time, tStep - tRec, tStep) / Unit::s;
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
	return result.reshape({size, 14UL});
}

double StarChi2Ptr(double *x, int metric, int gr_switch, int fix_switch, const py::array_t<double> &obs_time, const py::array_t<double> &obs_redshift, double redshift_sigma, const py::array_t<double> &obs_ra, double ra_sigma, const py::array_t<double> &obs_dec, double dec_sigma) {
	if (x[9] > 0.99)
		return GSL_POSINF;
	if (metric == 0 && gr_switch > 0) {
		x[7] = 1.;
		gr_switch = 0;
	}
	if (fix_switch & 1)
		x[1] = 0.07;
	if (fix_switch & 2)
		x[2] = -0.9;
	if (fix_switch & 4)
		x[4] = 0.0341;
	if (fix_switch & 8)
		x[5] = 0.080;
	if (fix_switch & 16)
		x[6] = -1.6;
	auto obs_data = CalculateStarOrbit(x[0], metric, x[7], x[3], x[13], x[8], x[9], x[10], x[11], x[12], gr_switch > 0, (gr_switch & 4) > 0, obs_time);
	double redshift_prob = 0., ra_prob = 0., dec_prob = 0.;
	const int size = obs_redshift.size();
	if (gr_switch & 1)
		for (int i = 0; i < size; ++i) {
			redshift_prob += gsl_pow_2(obs_redshift.at(i) - (obs_data.at(i, 12) - 1.) * 299792.458 - x[6]);
		}
	else
		for (int i = 0; i < size; ++i)
			redshift_prob += gsl_pow_2(obs_redshift.at(i) - obs_data.at(i, 9) - x[6]);
	if (gr_switch & 2)
		for (int i = 0; i < size; ++i) {
			ra_prob += gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 10) * x[0] * 9.870628713769018e-6 / x[3] - x[2] - (obs_time.at(i) - 2009.02) * x[5]);
			dec_prob += gsl_pow_2(obs_dec.at(i) - obs_data.at(i, 11) * x[0] * 9.870628713769018e-6 / x[3] - x[1] - (obs_time.at(i) - 2009.02) * x[4]);
		}
	else
		for (int i = 0; i < size; ++i) {
			ra_prob += gsl_pow_2(obs_ra.at(i) + obs_data.at(i, 2) * x[0] * 9.870628713769018e-6 / x[3] - x[2] - (obs_time.at(i) - 2009.02) * x[5]);
			dec_prob += gsl_pow_2(obs_dec.at(i) + obs_data.at(i, 1) * x[0] * 9.870628713769018e-6 / x[3] - x[1] - (obs_time.at(i) - 2009.02) * x[4]);
		}
	double chi2 = redshift_prob / gsl_pow_2(redshift_sigma) + ra_prob / gsl_pow_2(ra_sigma) + dec_prob / gsl_pow_2(dec_sigma);
	if (isnan(chi2))
		for (int i = 0; i < 14; ++i)
			py::print(x[i]);
	return chi2;
}

double StarChi2(py::array_t<double> x, int metric, int gr_switch, int fix_switch, py::array_t<double> obs_time, py::array_t<double> obs_redshift, double redshift_sigma, py::array_t<double> obs_ra, double ra_sigma, py::array_t<double> obs_dec, double dec_sigma) {
	return StarChi2Ptr(x.mutable_data(), metric, gr_switch, fix_switch, obs_time, obs_redshift, redshift_sigma, obs_ra, ra_sigma, obs_dec, dec_sigma);
}

py::array_t<double> GetStarChi2(py::array_t<double> mcmc_chain, int metric, int gr_switch, int fix_switch, py::array_t<double> obs_time, py::array_t<double> obs_redshift, double redshift_sigma, py::array_t<double> obs_ra, double ra_sigma, py::array_t<double> obs_dec, double dec_sigma) {
	const int m = mcmc_chain.shape(0);
	py::array_t<double> chi2(m);
	py::gil_scoped_release release;
#pragma omp parallel for
	for (int i = 0; i < m; ++i) {
		chi2.mutable_at(i) = StarChi2Ptr(mcmc_chain.mutable_data(i), metric, gr_switch, fix_switch, obs_time, obs_redshift, redshift_sigma, obs_ra, ra_sigma, obs_dec, dec_sigma);
	}
	py::gil_scoped_acquire acquire;
	return chi2;
}

py::array_t<double> HSExit(const py::array_t<double> &x) {
	for (int i = 0; i < x.size(); ++i)
		py::print(fmt::format("x[{}]={}", i, x.at(i)));
	return py::array_t<double>();
}

py::array_t<double> CalculateFullHSOrbit(const py::array_t<double> &x, int metric, int mode, bool ray_tracing, bool gr_time_delay, double t1) {
	Unit::Initialize(x.at(0)); // double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination,
	const double R = x.at(2) * Unit::pc, r = x.at(7) * R * Unit::mas, inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
	shared_ptr<Metric> main_metric;
	unique_ptr<View> view_ptr;
	if (metric == 1) {
		main_metric = make_shared<Schwarzschild>();
		view_ptr = make_unique<View>(make_unique<Schwarzschild>(), R, inclination, rotation);
	} else if (metric == 2) {
		main_metric = make_shared<Kerr>(x.at(1));
		view_ptr = make_unique<View>(make_shared<Kerr>(x.at(1)), R, inclination, rotation);
	} else {
		PrintlnError("Metric Error!");
		return py::array_t<double>();
	}
	time_system hotspot_time = T;
	HotSpot star_0(main_metric, hotspot_time, LAGRANGIAN, 0., 1., 1000., 1., false);
	size_t tStepNumber = 10000;
	double position[13], t = 0., tStep = 0., tRec = t1 / tStepNumber * Unit::s;
	if (mode == 0) { // circular
		if (star_0.InitializeCircular(r, x.at(8) * M_PI / 180., x.at(9)) != GSL_SUCCESS)
			return py::array_t<double>();
	} else if (mode == 1) { // geodesic
		if (star_0.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., x.at(11), x.at(12), inclination, rotation) != GSL_SUCCESS) {
			PrintlnError("Initialize Geodesic Error!");
			return HSExit(x);
		}
	} else if (mode == 2) { // helical
		if (star_0.InitializeHelical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11)) != GSL_SUCCESS) {
			PrintlnError("Initialize Helical Error!");
			return HSExit(x);
		}
	} else {
		PrintlnWarning("Mode Error!");
		return py::array_t<double>();
	}
	star_0.Position(position);
	if (gr_time_delay) {
		if (view_ptr->Trace(position, hotspot_time, position + 8, true) != GSL_SUCCESS) {
			PrintlnError("Initially trace star Error!");
			return HSExit(x);
		}
	}
	double h = 1.;
	int status = 0;
	auto result = py::array_t<double>(tStepNumber * 16);
	double *result_ptr = const_cast<double *>(result.data());
	for (size_t i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		if (status = star_0.IntegratorApply(&t, tStep, &h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		star_0.Position(result_ptr);
		if (gr_time_delay) {
			if (view_ptr->Trace(result_ptr, hotspot_time, result_ptr + 10, false) != GSL_SUCCESS) {
				PrintlnError("Trace star this position Error!");
				return HSExit(x);
			}
			result_ptr[15] = star_0.Luminosity(tStep, result_ptr[12]);
		}
		result_ptr[8] = t / Unit::s;
		const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
		result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
		result_ptr += 15;
	}
	return result.reshape({tStepNumber, 15UL});
}

py::array_t<double> CalculateHSOrbit(const py::array_t<double> &x, int metric, int mode, bool ray_tracing, bool gr_time_delay, const py::array_t<double> &obs_time) {
	const vector<double> estimate_step = {60., 30., 10., 3., 1., 0.5, 0.3, 0.2, 0.1, 0.};
	size_t estimate_idx = 0;
	double this_gr_obs_time, last_gr_obs_time, gr_offset = 0.;
	Unit::Initialize(x.at(0)); // double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination,
	double R = x.at(2) * Unit::pc, r = x.at(7) * R * Unit::mas, inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
	shared_ptr<Metric> main_metric;
	unique_ptr<View> view_ptr;
	if (metric == 1) {
		main_metric = make_shared<Schwarzschild>();
		view_ptr = make_unique<View>(make_unique<Schwarzschild>(), R, inclination, rotation);
	} else if (metric == 2) {
		main_metric = make_shared<Kerr>(x.at(1));
		view_ptr = make_unique<View>(make_shared<Kerr>(x.at(1)), R, inclination, rotation);
	} else {
		PrintlnError("Metric Error!");
		return py::array_t<double>();
	}
	time_system hotspot_time = T;
	HotSpot star_0(main_metric, hotspot_time, LAGRANGIAN, 0., 1., 1000., 1., false);
	double *last_position = new double[14], *this_position = new double[14], t = 0., tStep = 0., tRec = 0.1 * Unit::s, t0, last_obs_time = 0., this_obs_time;
	const double sin_inc = sin(inclination), cos_inc = cos(inclination);
	if (mode == 0) { // circular
		if (star_0.InitializeCircular(r, x.at(8) * M_PI / 180., x.at(9)) != GSL_SUCCESS)
			return py::array_t<double>();
	} else if (mode == 1) { // geodesic
		if (star_0.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., x.at(11), x.at(12), inclination, rotation) != GSL_SUCCESS)
			return py::array_t<double>();
	} else if (mode == 2) { // helical
		if (star_0.InitializeHelical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11)) != GSL_SUCCESS)
			return py::array_t<double>();
	} else {
		PrintlnWarning("Mode Error!");
		return py::array_t<double>();
	}
	star_0.Position(last_position);
	double x0 = last_position[1] * abs(sin(last_position[2])) * cos(last_position[3]), z0 = last_position[1] * GSL_SIGN(last_position[2]) * cos(last_position[2]);
	if (gr_time_delay) {
		if (view_ptr->Trace(last_position, hotspot_time, last_position + 8, false) != GSL_SUCCESS) {
			PrintlnError("Initially trace star Error!");
			return HSExit(x);
		}
		t0 = last_position[11];
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
	auto result = py::array_t<double>(size * 16);
	double *result_ptr = const_cast<double *>(result.data());
	for (int i = 0;; ++i) {
		tStep += tRec;
		if (status = star_0.IntegratorApply(&t, tStep, &h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		star_0.Position(this_position);
		if (this_position[1] < 4.)
			return py::array_t<double>();
		this_obs_time = (t + vz0 - this_position[1] * (abs(sin(this_position[2])) * cos(this_position[3]) * sin_inc + GSL_SIGN(this_position[2]) * cos(this_position[2]) * cos_inc)) / Unit::s;
		if (gr_time_delay) {
			if (this_obs_time + gr_offset + estimate_step[estimate_idx] > obs_time.at(idx)) {
				if (++estimate_idx < estimate_step.size()) {
					if (view_ptr->Trace(last_position, hotspot_time, last_position + 8, true) != GSL_SUCCESS) {
						PrintlnError("Trace star last position Error!");
						return HSExit(x);
					}
					last_gr_obs_time = t / Unit::s + t0 - last_position[11];
					last_position[13] = star_0.Luminosity(tStep, last_position[10]);
					gr_offset = last_gr_obs_time - this_obs_time;
				} else {
					estimate_idx = 0;
					if (view_ptr->Trace(this_position, hotspot_time, this_position + 8, true) != GSL_SUCCESS) {
						PrintlnError("Trace star this position Error!");
						return HSExit(x);
					}
					this_gr_obs_time = t / Unit::s + t0 - this_position[11];
					this_position[13] = star_0.Luminosity(tStep, this_position[10]);
					gr_offset = 0.;
					InterpolateSphericalPositionToCartesian(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, last_position, this_position, result_ptr);
					result_ptr[8] = LinearInterpolation(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, tStep - tRec, tStep) / Unit::s;
					LinearInterpolation(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, last_position + 8, this_position + 8, result_ptr + 10, 6);
					const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
					result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
					result_ptr += 15;
					if (++idx >= size)
						break;
				}
			}
		} else if (this_obs_time > obs_time.at(idx)) {
			InterpolateSphericalPositionToCartesian(obs_time.at(idx), last_obs_time, this_obs_time, last_position, this_position, result_ptr);
			if (ray_tracing || metric == 2) {
				CartesianToSpherical(result_ptr);
				if (ray_tracing) {
					if (view_ptr->Trace(result_ptr, hotspot_time, result_ptr + 10, true) != GSL_SUCCESS) {
						PrintlnError("Trace star position Error!");
						return HSExit(x);
					}
					result_ptr[15] = star_0.Luminosity(tStep, result_ptr[12]);
				}
				if (metric == 2)
					result_ptr[1] -= 1.;
				SphericalToCartesian(result_ptr);
			}
			result_ptr[8] = LinearInterpolation(obs_time.at(idx), last_obs_time, this_obs_time, tStep - tRec, tStep) / Unit::s;
			const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
			result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
			result_ptr += 15;
			if (++idx >= size)
				break;
		}
		last_obs_time = this_obs_time;
		swap(last_position, this_position);
	}
	return result.reshape({size, 15UL});
}

double HSChi2(py::array_t<double> &x, int metric, int mode, int gr_switch, int fix_switch, py::array_t<double> &obs_time, py::array_t<double> &obs_flux, py::array_t<double> &flux_sigma, py::array_t<double> &obs_ra, py::array_t<double> &ra_sigma, py::array_t<double> &obs_dec, py::array_t<double> &dec_sigma) {
	if (metric == 0 && gr_switch > 0) {
		x.mutable_at(1) = 1.;
		gr_switch = 0;
	}
	if (mode == 1 && fix_switch & 4096) {
		Unit::Initialize(x.at(0));
		const double g00 = 1. - 2. / (x.at(7) * x.at(2) * Unit::pc * Unit::mas);
		// dt/dtau = 1./g00;
		// y[4] == g00
		// g00 - g00 * g00 = y[5] * y[5] / g00 + v_phi * v_phi;
		x.mutable_at(12) = sqrt(g00 - gsl_pow_2(g00) - gsl_pow_2(x.at(11)) / g00);
	}
	auto obs_data = CalculateHSOrbit(x, metric, mode, gr_switch > 0, (gr_switch & 4) > 0, obs_time);
	if (obs_data.size() == 0)
		return GSL_POSINF;
	double flux_prob = 0., ra_prob = 0., dec_prob = 0.;
	const int size = obs_flux.size();
	if (mode == 0 && !(fix_switch & 1024))
		for (int i = 0; i < size; ++i)
			flux_prob += gsl_pow_2((obs_flux.at(i) - abs(obs_data.at(i, 14) * obs_data.at(i, 15))) / flux_sigma.at(i));
	if (gr_switch & 2)
		for (int i = 0; i < size; ++i) {
			ra_prob += gsl_pow_2((obs_ra.at(i) + obs_data.at(i, 10) * x.at(0) * 9.870628713769018e-6 / x.at(2) - x.at(3)) / ra_sigma.at(i));
			dec_prob += gsl_pow_2((obs_dec.at(i) - obs_data.at(i, 11) * x.at(0) * 9.870628713769018e-6 / x.at(2) - x.at(4)) / dec_sigma.at(i));
		}
	else {
		PrintlnError("Inclination and rotation is ignored!");
		for (int i = 0; i < size; ++i) {
			ra_prob += gsl_pow_2((obs_ra.at(i) + obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(2) - x.at(3)) / ra_sigma.at(i));	   // FIXME: inclination!
			dec_prob += gsl_pow_2((obs_dec.at(i) + obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(2) - x.at(4)) / dec_sigma.at(i)); // FIXME: inclination!
		}
	}
	return flux_prob + ra_prob + dec_prob;
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

	m.def("get_max_threads", &omp_get_max_threads, "Returns max number of threads");
	m.def("set_num_threads", &omp_set_num_threads, "Set number of threads");

	m.def("CalculatePericenterTime", &CalculatePericenterTime, R"pbdoc(
		Get observed pericenter time.
	)pbdoc");

	m.def("CalculateFullStarOrbit", &CalculateFullStarOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("CalculateStarOrbit", &CalculateStarOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("StarChi2", &StarChi2, R"pbdoc(
		Get chi2
	)pbdoc");

	m.def("GetStarChi2", &GetStarChi2, R"pbdoc(
		Get chi2
	)pbdoc");

	m.def("CalculateFullHSOrbit", &CalculateFullHSOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("CalculateHSOrbit", &CalculateHSOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("HSChi2", &HSChi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

	m.def("EllipticIntegral", &EllipticIntegral, R"pbdoc(
		Elliptic Integral
	)pbdoc");

	m.def("EllipticIntegral2Complex", &EllipticIntegral2Complex, R"pbdoc(
		Elliptic Integral with 2 complex roots
	)pbdoc");

	m.def("EllipticIntegral4Complex", &EllipticIntegral4Complex, R"pbdoc(
		Elliptic Integral with 4 complex roots
	)pbdoc");

	m.attr("s") = Unit::s;
	m.attr("day") = Unit::day;
	m.attr("yr") = Unit::yr;
	m.attr("__version__") = VERSION;
}
