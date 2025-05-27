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

constexpr auto HOTSPOT_RETURN_WIDTH = 17UL;
double CalculatePericenterTime(double mass, int metric, double fSP, double R, double t_apo, double a, double e, double inclination, double ascending_node, double periapsis) {
#ifdef GSL_RANGE_CHECK_OFF
	gsl_set_error_handler_off();
#endif
	Unit::Initialize(mass);
	R *= Unit::pc;
	a *= R * Unit::mas;
	inclination *= M_PI / 180.;
	ascending_node *= M_PI / 180.;
	periapsis *= M_PI / 180.;
	double t = (t_apo - 2002.) * Unit::yr, tStep = 16.3 * Unit::yr, tRec = 0.000001 * Unit::yr;
	shared_ptr<Metric> main_metric;
	unique_ptr<View> view_ptr;
	if (metric == 0)
		main_metric = make_shared<PN1>(fSP);
	else {
		main_metric = make_shared<Schwarzschild>();
		view_ptr = make_unique<View>(make_unique<Schwarzschild>(), R, 0., 0.);
	}
	Particle star_0(main_metric, T, LAGRANGIAN, false);
	if (metric == 2)
		star_0.InitializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
	else
		star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
	double h = 1., position[8], last_radius = GSL_POSINF, t0, record[4];
	star_0.Position(position);
	if (metric == 0)
		t0 = position[1] * GSL_SIGN(position[2]) * cos(position[2]);
	else {
		view_ptr->Trace(position, T, record, false);
		t0 = record[3];
	}
	int status = 0;
	for (int i = 0;; ++i) {
		status = star_0.IntegratorApply(&t, tStep, &h);
		if (status > 0)
			PrintlnError("main status = {}", status);
		star_0.Position(position);
		if (position[1] > last_radius) {
			if (metric == 0)
				return (tStep + t0 - position[1] * GSL_SIGN(position[2]) * cos(position[2])) / Unit::yr + 2002.;
			else {
				view_ptr->Trace(position, T, record, false);
				return (tStep + (t0 - record[3]) * Unit::s) / Unit::yr + 2002.;
			}
		}
		last_radius = position[1];
		tStep += tRec;
	}
}

py::array_t<double> CalculateFullStarOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e, double inclination, double ascending_node, double periapsis, double t1) {
#ifdef GSL_RANGE_CHECK_OFF
	gsl_set_error_handler_off();
#endif
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
	} else if (metric == 3) {
		main_metric = make_shared<Newton>(3);
	} else {
		main_metric = make_shared<Schwarzschild>();
		view_ptr = make_unique<View>(make_unique<Schwarzschild>(), R, 0., 0.);
	}
	TimeSystem star_time = T;
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
#ifdef GSL_RANGE_CHECK_OFF
	gsl_set_error_handler_off();
#endif
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
	} else if (metric == 3) {
		main_metric = make_shared<Newton>(3);
	} else {
		main_metric = make_shared<Schwarzschild>();
		view_ptr = make_unique<View>(make_unique<Schwarzschild>(), R, M_PI_4, 0.);
	}
	TimeSystem star_time = T;
	Particle star_0(main_metric, star_time, LAGRANGIAN, false);
	double z0, t0, last_obs_time = 2002., this_obs_time, *last_position = new double[13], *this_position = new double[13];
	if (metric == 2)
		star_0.InitializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, 0., M_PI_4, 0.);
	else
		star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, 0., M_PI_4, 0.);
	star_0.Position(last_position);
	z0 = last_position[1] * GSL_SIGN(last_position[2]) * cos(last_position[2]);
	if (gr_time_delay) {
		if (view_ptr->Trace(last_position, star_time, last_position + 8, false) != GSL_SUCCESS) {
			PrintlnError("Initially trace star Error!");
			return py::array_t<double>();
		}
		t0 = last_position[11];
	}
	if (metric == 2)
		star_0.InitializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, M_PI, M_PI_4, 0.);
	else
		star_0.InitializeKeplerian(a, e, inclination, periapsis, ascending_node, M_PI, M_PI_4, 0.);
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
				if (view_ptr->Trace(last_position, star_time, last_position + 8, false) != GSL_SUCCESS) {
					PrintlnError("Trace last_position Error!");
					return py::array_t<double>();
				}
				if (view_ptr->Trace(this_position, star_time, this_position + 8, false) != GSL_SUCCESS) {
					PrintlnError("Trace this_position Error!");
					return py::array_t<double>();
				}
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
						if (view_ptr->Trace(result_ptr, star_time, result_ptr + 10, false) != GSL_SUCCESS) {
							PrintlnError("Trace star position Error!");
							return py::array_t<double>();
						}
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

double StarChi2(py::array_t<double> x, int metric, int gr_switch, py::array_t<double> obs_time, py::array_t<double> obs_redshift, py::array_t<double> redshift_sigma, py::array_t<double> obs_ra, py::array_t<double> ra_sigma, py::array_t<double> obs_dec, py::array_t<double> dec_sigma) {
	if (x.at(9) > 0.99)
		return GSL_POSINF;
	if (metric == 0 && gr_switch > 0) {
		x.mutable_at(7) = 1.;
		gr_switch = 0;
	}
	auto obs_data = CalculateStarOrbit(x.at(0), metric, x.at(7), x.at(3), x.at(13), x.at(8), x.at(9), x.at(10), x.at(11), x.at(12), gr_switch > 0, (gr_switch & 4) > 0, obs_time);
	double redshift_prob = 0., ra_prob = 0., dec_prob = 0.;
	const int size = obs_redshift.size();
	if (gr_switch & 1) {
		for (int i = 0; i < size; ++i)
			if (obs_redshift.at(i) != GSL_POSINF)
				redshift_prob += gsl_pow_2((obs_redshift.at(i) - (obs_data.at(i, 12) - 1.) * 299792.458 - x.at(6)) / redshift_sigma.at(i));
	} else
		for (int i = 0; i < size; ++i)
			if (obs_redshift.at(i) != GSL_POSINF)
				redshift_prob += gsl_pow_2((obs_redshift.at(i) - obs_data.at(i, 9) - x.at(6)) / redshift_sigma.at(i));
	if (gr_switch & 2)
		for (int i = 0; i < size; ++i) {
			if (obs_ra.at(i) != GSL_POSINF)
				ra_prob += gsl_pow_2((obs_ra.at(i) + obs_data.at(i, 10) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - (obs_time.at(i) - 2009.02) * x.at(5)) / ra_sigma.at(i));
			if (obs_dec.at(i) != GSL_POSINF)
				dec_prob += gsl_pow_2((obs_dec.at(i) - obs_data.at(i, 11) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - (obs_time.at(i) - 2009.02) * x.at(4)) / dec_sigma.at(i));
		}
	else
		for (int i = 0; i < size; ++i) {
			if (obs_ra.at(i) != GSL_POSINF)
				ra_prob += gsl_pow_2((obs_ra.at(i) + obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(2) - (obs_time.at(i) - 2009.02) * x.at(5)) / ra_sigma.at(i));
			if (obs_dec.at(i) != GSL_POSINF)
				dec_prob += gsl_pow_2((obs_dec.at(i) + obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(3) - x.at(1) - (obs_time.at(i) - 2009.02) * x.at(4)) / dec_sigma.at(i));
		}
	double chi2 = redshift_prob + ra_prob + dec_prob;
	if (isnan(chi2)) {
		if (isnan(redshift_prob))
			py::print("[!] redshift_prob is NaN");
		if (isnan(ra_prob))
			py::print("[!] ra_prob is NaN");
		if (isnan(dec_prob))
			py::print("[!] dec_prob is NaN");
		auto obs_data_test = CalculateStarOrbit(x.at(0), metric, x.at(7), x.at(3), x.at(13), x.at(8), x.at(9), x.at(10), x.at(11), x.at(12), gr_switch > 0, (gr_switch & 4) > 0, obs_time);
		for (int i = 0; i < 14; ++i)
			py::print(x.at(i));
	}
	return chi2;
}

py::array_t<double> HSExit(const py::array_t<double> &x) {
	return py::array_t<double>();
	for (int i = 0; i < x.size(); ++i)
		py::print(fmt::format("x[{}]={}", i, x.at(i)));
	return py::array_t<double>();
}

py::array_t<double> CalculateFullHSOrbit(const py::array_t<double> &x, int metric, int mode, bool ray_tracing, bool gr_time_delay, double t1) {
#ifdef GSL_RANGE_CHECK_OFF
	gsl_set_error_handler_off();
#endif
	Unit::Initialize(x.at(0)); // double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination,
	const double R = x.at(2) * Unit::pc, r = x.at(7) * R * Unit::mas, inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
	const int offset[] = {10, 13, 12};
	shared_ptr<Metric> main_metric;
	unique_ptr<View> view_ptr;
	if (metric == 1 || (metric == 2 && abs(x.at(1)) < GSL_ROOT4_DBL_EPSILON)) {
		main_metric = make_shared<Schwarzschild>();
	} else if (metric == 2) {
		main_metric = make_shared<Kerr>(x.at(1));
	} else {
		PrintlnError("Metric Error!");
		return py::array_t<double>();
	}
	view_ptr = make_unique<View>(main_metric, R, inclination, rotation);
	TimeSystem hotspot_time = T;
	HotSpot hotspot(main_metric, hotspot_time, LAGRANGIAN, x.at(offset[mode]), x.at(offset[mode] + 1), x.at(offset[mode] + 2) * Unit::s, x.at(offset[mode] + 3), false);
	size_t tStepNumber = 200;
	double position[13], t = 0., tStep = 0., tRec = t1 / tStepNumber * Unit::s;
	if (mode == 0) { // circular
		if (hotspot.InitializeCircular(r, x.at(8) * M_PI / 180., x.at(9)) != GSL_SUCCESS) {
			PrintlnError("Initialize Circular Error!");
			return HSExit(x);
		}
	} else if (mode == 1) { // geodesic
		if (hotspot.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., x.at(11), x.at(12), inclination, rotation) != GSL_SUCCESS) {
			PrintlnError("Initialize Geodesic Error!");
			return HSExit(x);
		}
	} else if (mode == 2) { // helical
		if (hotspot.InitializeHelical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11)) != GSL_SUCCESS) {
			PrintlnError("Initialize Helical Error!");
			return HSExit(x);
		}
	} else {
		PrintlnWarning("Mode Error!");
		return py::array_t<double>();
	}
	hotspot.Position(position);
	if (gr_time_delay) {
		if (view_ptr->Trace(position, hotspot_time, position + 8, true) != GSL_SUCCESS) {
			PrintlnError("Initially trace star Error!");
			return HSExit(x);
		}
	}
	double h = 1.;
	int status = 0;
	auto result = py::array_t<double>(tStepNumber * HOTSPOT_RETURN_WIDTH);
	double *result_ptr = const_cast<double *>(result.data());
	for (size_t i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		if (status = hotspot.IntegratorApply(&t, tStep, &h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		hotspot.Position(result_ptr);
		if (gr_time_delay) {
			if (view_ptr->Trace(result_ptr, hotspot_time, result_ptr + 10, true) != GSL_SUCCESS) {
				PrintlnError("Trace star this position Error!");
				return HSExit(x);
			}
			result_ptr[15] = Flux(hotspot.Luminosity(tStep), result_ptr[14], result_ptr[12]);
			result_ptr[16] = FluxDensity(hotspot.SpectralDensity(tStep, result_ptr[12]), result_ptr[14]);
		}
		SphericalToCartesian(result_ptr);
		result_ptr[8] = t / Unit::s;
		const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
		result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
		result_ptr += HOTSPOT_RETURN_WIDTH;
	}
	return result.reshape({tStepNumber, HOTSPOT_RETURN_WIDTH});
}

py::array_t<double> CalculateHSOrbit(const py::array_t<double> &x, int metric, int mode, bool ray_tracing, bool gr_time_delay, const py::array_t<double> &obs_time, bool calculate_magnification) {
#ifdef GSL_RANGE_CHECK_OFF
	gsl_set_error_handler_off();
#endif
	auto t_start = chrono::steady_clock::now();
	const vector<double> estimate_step = {60., 30., 10., 3., 1., 0.5, 0.1, 0.};
	size_t estimate_idx = 0;
	double this_gr_obs_time, last_gr_obs_time, gr_offset = 0.;
	Unit::Initialize(x.at(0)); // double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination,
	double R = x.at(2) * Unit::pc, r = x.at(7) * R * Unit::mas, inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
	const int offset[] = {10, 13, 12};
	shared_ptr<Metric> main_metric;
	unique_ptr<View> view_ptr;
	if (metric == 1 || (metric == 2 && abs(x.at(1)) < GSL_ROOT4_DBL_EPSILON))
		main_metric = make_shared<Schwarzschild>();
	else if (metric == 2)
		main_metric = make_shared<Kerr>(x.at(1));
	else {
		PrintlnError("Metric Error!");
		return py::array_t<double>();
	}
	view_ptr = make_unique<View>(main_metric, R, inclination, rotation);
	TimeSystem hotspot_time = T;
	HotSpot hotspot(main_metric, hotspot_time, LAGRANGIAN, x.at(offset[mode]), x.at(offset[mode] + 1), x.at(offset[mode] + 2) * Unit::s, x.at(offset[mode] + 3) * Unit::s, false);
	double t = 0., tStep = 0., tRec = 0.1 * Unit::s, t0, last_obs_time = 0., this_obs_time;
	array<double, 8> last_position, this_position;
	array<double, 7> last_record, this_record;
	const double sin_inc = sin(inclination), cos_inc = cos(inclination);
	if (mode == 0) { // circular
		if (hotspot.InitializeCircular(r, x.at(8) * M_PI / 180., x.at(9)) != GSL_SUCCESS)
			return HSExit(x);
	} else if (mode == 1) { // geodesic
		if (hotspot.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., x.at(11), x.at(12), inclination, rotation) != GSL_SUCCESS)
			return HSExit(x);
	} else if (mode == 2) { // helical
		if (hotspot.InitializeHelical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11)) != GSL_SUCCESS)
			return HSExit(x);
	} else {
		PrintlnWarning("Mode Error!");
		return py::array_t<double>();
	}
	hotspot.Position(last_position.data());
	if (last_position[1] < 3.) {
		// py::print("[!] R < 3.0, stop!");
		return py::array_t<double>();
	}
	double x0 = last_position[1] * abs(sin(last_position[2])) * cos(last_position[3]), z0 = last_position[1] * GSL_SIGN(last_position[2]) * cos(last_position[2]);
	if (gr_time_delay) {
		if (view_ptr->Trace(last_position.data(), hotspot_time, last_record.data(), calculate_magnification) != GSL_SUCCESS) {
			PrintlnError("Initially trace star Error!");
			return HSExit(x);
		}
		t0 = last_record[3];
	}
	double vz0 = x0 * sin_inc + z0 * cos_inc;
	double h = -1.;
	int status = 0;
	if (status = hotspot.IntegratorApply(&t, 0., &h); status != 0)
		py::print("[!] IntegratorApply status =", status);
	h = 1.;
	if (status = hotspot.IntegratorReset(); status != 0)
		py::print("[!] IntegratorReset status =", status);
	size_t idx = 0, size = obs_time.size();
	auto result = py::array_t<double>(size * HOTSPOT_RETURN_WIDTH);
	double *result_ptr = const_cast<double *>(result.data());
	for (int i = 0;; ++i) {
#ifdef GSL_RANGE_CHECK_OFF
		if (auto now = chrono::steady_clock::now(); now - t_start > chrono::seconds(10)) {
			py::print("CalculateHSOrbit() timeout!");
			return HSExit(x);
		}
#endif
		tStep += tRec;
		if (status = hotspot.IntegratorApply(&t, tStep, &h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		hotspot.Position(this_position.data());
		if (this_position[1] < 3.) {
			// py::print("[!] R < 3.0, stop!");
			return py::array_t<double>();
		}
		this_obs_time = (t + vz0 - this_position[1] * (abs(sin(this_position[2])) * cos(this_position[3]) * sin_inc + GSL_SIGN(this_position[2]) * cos(this_position[2]) * cos_inc)) / Unit::s;
		if (gr_time_delay) {
			if (this_obs_time + gr_offset + estimate_step[estimate_idx] > obs_time.at(idx)) {
				if (view_ptr->Trace(this_position.data(), hotspot_time, this_record.data(), calculate_magnification) != GSL_SUCCESS) {
					if (this_obs_time + gr_offset - 0.1 > obs_time.at(idx)) {
						PrintlnError("Trace star this position Error!");
						return HSExit(x);
					}
					if (estimate_idx + 2 < estimate_step.size())
						++estimate_idx;
					continue;
				}
				this_gr_obs_time = t / Unit::s + t0 - this_record[3];
				if (calculate_magnification) {
					this_record[5] = Flux(hotspot.Luminosity(tStep), this_record[4], this_record[2]);
					this_record[6] = FluxDensity(hotspot.SpectralDensity(tStep, this_record[2]), this_record[4]);
				}
				if (gr_offset = this_gr_obs_time - this_obs_time; abs(gr_offset) > 3600.) { // difference between Romer delay and full GR delay
					PrintlnError("gr_offset > 3600 s!");
					for (int i = 0; i < 8; ++i)
						py::print("this_position[", i, "]=", this_position[i]);
					return HSExit(x);
				}
				if (++estimate_idx >= estimate_step.size()) {
					estimate_idx = 0;
					InterpolateSphericalPositionToCartesian(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, last_position.data(), this_position.data(), result_ptr);
					result_ptr[8] = LinearInterpolation(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, tStep - tRec, tStep) / Unit::s;
					const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
					result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
					LinearInterpolation(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, last_record.data(), this_record.data(), result_ptr + 10, 7);
					result_ptr += HOTSPOT_RETURN_WIDTH;
					if (++idx >= size)
						break;
				}
				last_gr_obs_time = this_gr_obs_time;
				swap(last_record, this_record);
			}
		} else if (this_obs_time > obs_time.at(idx)) {
			InterpolateSphericalPositionToCartesian(obs_time.at(idx), last_obs_time, this_obs_time, last_position.data(), this_position.data(), result_ptr);
			if (ray_tracing) {
				CartesianToSpherical(result_ptr);
				if (ray_tracing) {
					if (view_ptr->Trace(result_ptr, hotspot_time, result_ptr + 10, calculate_magnification) != GSL_SUCCESS) {
						PrintlnError("Trace star position Error!");
						return HSExit(x);
					}
					if (calculate_magnification) {
						result_ptr[15] = Flux(hotspot.Luminosity(tStep), result_ptr[14], result_ptr[12]);
						result_ptr[16] = FluxDensity(hotspot.SpectralDensity(tStep, result_ptr[12]), result_ptr[14]);
					}
				}
				SphericalToCartesian(result_ptr);
			}
			result_ptr[8] = LinearInterpolation(obs_time.at(idx), last_obs_time, this_obs_time, tStep - tRec, tStep) / Unit::s;
			const double delta_epsilon = 1. - 2. / Norm(result_ptr + 1);
			result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - Dot(result_ptr + 5)) - 1.) * 299792.458;
			result_ptr += HOTSPOT_RETURN_WIDTH;
			if (++idx >= size)
				break;
		}
		last_obs_time = this_obs_time;
		swap(last_position, this_position);
	}
	return result.reshape({size, HOTSPOT_RETURN_WIDTH});
}

double HSFluxChi2(py::array_t<double> &x, py::array_t<double> &intrinsic_time, py::array_t<double> &obs_magnification, py::array_t<double> &obs_redshift, py::array_t<double> &obs_flux, py::array_t<double> &flux_sigma) {
	const int size = obs_flux.size();
	shared_ptr<Metric> main_metric;
	HotSpot hotspot(main_metric, TAU, LAGRANGIAN, x.at(0), x.at(1), x.at(2) * Unit::s, x.at(3) * Unit::s, false);
	double chi2 = 0.;
	for (int i = 0; i < size; ++i)
		chi2 += gsl_pow_2((FluxDensity(hotspot.SpectralDensity(intrinsic_time.at(i), obs_redshift.at(i)), obs_magnification.at(i)) - obs_flux.at(i)) / flux_sigma.at(i));
	return chi2;
}
py::array_t<double> HSDistinctChi2(py::array_t<double> &x, int metric, int mode, int gr_switch, py::array_t<double> &obs_time, py::array_t<double> &obs_flux, py::array_t<double> &flux_sigma, py::array_t<double> &obs_ra, py::array_t<double> &ra_sigma, py::array_t<double> &obs_dec, py::array_t<double> &dec_sigma, double disk_flux) {
	if (metric == 0 && gr_switch > 0) {
		x.mutable_at(1) = 1.;
		gr_switch = 0;
	}
	auto obs_data = CalculateHSOrbit(x, metric, mode, gr_switch > 0, (gr_switch & 4) > 0, obs_time, gr_switch & 8);
	if (obs_data.size() == 0)
		return py::array_t<double>();
	const int size = obs_flux.size();
	auto result = py::array_t<double>(size);
	double *result_ptr = const_cast<double *>(result.data());
	memset(result_ptr, 0, size * sizeof(double));
	if (gr_switch & 8)
		for (int i = 0; i < size; ++i)
			result_ptr[i] += gsl_pow_2((obs_flux.at(i) - obs_data.at(i, 16)) / flux_sigma.at(i));
	if (gr_switch & 2)
		for (int i = 0; i < size; ++i) { // 9.870628713769018e-6 = 1. / (Unit::pc * Unit::mas)
			result_ptr[i] += gsl_pow_2((obs_ra.at(i) + obs_data.at(i, 10) * x.at(0) * 9.870628713769018e-6 / x.at(2) * (1. - disk_flux / obs_flux.at(i)) - x.at(3)) / ra_sigma.at(i));
			result_ptr[i] += gsl_pow_2((obs_dec.at(i) - obs_data.at(i, 11) * x.at(0) * 9.870628713769018e-6 / x.at(2) * (1. - disk_flux / obs_flux.at(i)) - x.at(4)) / dec_sigma.at(i));
		}
	else {
		PrintlnError("Inclination and rotation is ignored!");
		for (int i = 0; i < size; ++i) {
			result_ptr[i] += gsl_pow_2((obs_ra.at(i) + obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(2) * (1. - disk_flux / obs_flux.at(i)) - x.at(3)) / ra_sigma.at(i));	// FIXME: inclination!
			result_ptr[i] += gsl_pow_2((obs_dec.at(i) + obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(2) * (1. - disk_flux / obs_flux.at(i)) - x.at(4)) / dec_sigma.at(i)); // FIXME: inclination!
		}
	}
	return result;
}
double HSChi2(py::array_t<double> &x, int metric, int mode, int gr_switch, py::array_t<double> &obs_time, py::array_t<double> &obs_flux, py::array_t<double> &flux_sigma, py::array_t<double> &obs_ra, py::array_t<double> &ra_sigma, py::array_t<double> &obs_dec, py::array_t<double> &dec_sigma) {
	cerr.rdbuf(cout.rdbuf());
	auto result = HSDistinctChi2(x, metric, mode, gr_switch, obs_time, obs_flux, flux_sigma, obs_ra, ra_sigma, obs_dec, dec_sigma, 0.0);
	if (result.size() == 0)
		return GSL_POSINF;
	double prob_sum = 0.0;
	for (int i = 0; i < result.size(); ++i)
		prob_sum += result.at(i);
	if (gr_switch & 16)
		prob_sum += (gsl_pow_2(x.at(3)) + gsl_pow_2(x.at(4))) * 400.;
	if (isnan(prob_sum)) {
		for (int i = 0; i < result.size(); ++i)
			py::print(result.at(i));
		return GSL_POSINF;
	}
	return prob_sum;
}

py::array_t<double> FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_target, const double theta_target, const double phi_target) {
	auto photon = py::array_t<double>(9);
	double *photon_ptr = const_cast<double *>(photon.data());
	double alpha = 0., beta = 0.;
	Schwarzschild().FastTrace(r_observer, theta_observer, sin_theta_observer, cos_theta_observer, r_target, theta_target, phi_target, alpha, beta, photon_ptr);
	return photon;
}

py::array_t<double> CppSort(py::array_t<double> data) {
	double *data_ptr = const_cast<double *>(data.data());
	sort(data_ptr, data_ptr + data.size());
	return data;
}

py::array_t<double> KerrCoverage(double spin, double theta) {
	// double u_roots[4];
	// PolySolveQuartic(-2225768.312261119, -9.4005613541838091 * -115123.80934897072, 0., -115123.80934897072, u_roots);
	Unit::Initialize(4.15e6);
	auto main_metric = Kerr(spin);
	const double sin_theta = sin(theta), cos_theta = cos(theta);
	vector<double> u(201);
	auto result = py::array_t<double>(201 * 201 * 4 * 201);
	memset(result.mutable_data(), 0xff, sizeof(double) * 201 * 201 * 4 * 201);
	if (!isnan(result.data()[0]))
		return py::array_t<double>();
	for (int i = 0; i < 201; ++i)
		u[i] = 1. / (7.3090991130925094 - 0.01 * i);
	int offset = 4 * 201;
	int offset_0 = 0;
	int offset_1 = 201;
	int offset_2 = offset_1 * 2;
	int offset_3 = offset_1 * 3;
	for (int i = -100; i < 101; ++i) {
		double alpha = 5.630656 + 0.0001 * i;
		for (int j = -100; j < 101; ++j) {
			double beta = 6.527776 + 0.0001 * j;
			main_metric.CalcThetaPhi(8.18 * Unit::pc, theta, sin_theta, cos_theta, alpha, beta, u, result.mutable_data() + offset_0, result.mutable_data() + offset_1, result.mutable_data() + offset_2, result.mutable_data() + offset_3);
			offset_0 += offset;
			offset_1 += offset;
			offset_2 += offset;
			offset_3 += offset;
			ProgressBar::bars_[0].set_progress((i + 100 + (j + 100) / 201.) * 100. / 201.);
		}
	}
	return result.reshape({201, 201, 4, 201});
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

	m.def("CalculateFullHSOrbit", &CalculateFullHSOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("CalculateHSOrbit", &CalculateHSOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("HSChi2", &HSChi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

	m.def("HSDistinctChi2", &HSDistinctChi2, R"pbdoc(
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

	m.def("FastTrace", &FastTrace, R"pbdoc(
		Schwarzschild fast trace
	)pbdoc");

	m.def("CppSort", &CppSort, R"pbdoc(
		C++ sort
	)pbdoc");

	m.def("KerrCoverage", &KerrCoverage, R"pbdoc(
		KerrCoverage
	)pbdoc");

	m.attr("s") = Unit::s;
	m.attr("day") = Unit::day;
	m.attr("yr") = Unit::yr;
	m.attr("__version__") = VERSION;
}
