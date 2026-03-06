/**
 * @file hotspot.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <array>
#include <cassert>
#include <chrono>
#include <csignal>
#include <filesystem>
#include <memory>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <boost/algorithm/algorithm.hpp>
#include <fmt/core.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include "infra/status.hpp"
#include "math/geometry.hpp"
#include "math/interpolation.hpp"
#include "metric/kerr.hpp"
#include "metric/schwarzschild.hpp"
#include "object/object.hpp"
#include "object/particle.hpp"
#include "units.hpp"
#include "view/flux.hpp"
#include "view/view.hpp"

namespace py = pybind11;
namespace ublas = boost::numeric::ublas;
using namespace std;
using namespace sbody;

#define PRECISION double

constexpr auto HOTSPOT_RETURN_WIDTH = 17UL;

py::array_t<double> HSExit(const py::array_t<double> &x) {
	return py::array_t<double>();
	for (int i = 0; i < x.size(); ++i)
		py::print(fmt::format("x[{}]={}", i, x.at(i)));
	return py::array_t<double>();
}

py::array_t<double> CalculateFullHSOrbit(const py::array_t<double> &x, int metric, int mode, bool ray_tracing, bool gr_time_delay, double t1) {
	Unit unit(x.at(0)); // double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination,
	const double R = x.at(2) * unit.pc, r = x.at(7) * R * unit.mas, inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
	const array<int, 8> offset = {10, 13, 12, 12, 12, 12, 13, 12};
	unique_ptr<Metric<double>> main_metric;
	unique_ptr<View<double>> view_ptr;
	if (metric == 1 || (metric == 2 && abs(x.at(1)) < boost::math::tools::forth_root_epsilon<double>())) {
		main_metric = make_unique<Schwarzschild<double>>();
	} else if (metric == 2) {
		main_metric = make_unique<Kerr<double>>(x.at(1));
	} else {
		PrintlnError("Metric Error!");
		return py::array_t<double>();
	}
	view_ptr = make_unique<View<double>>(*main_metric, R, inclination, rotation);
	TimeSystem hotspot_time = TimeSystem::T;
	Hotspot hotspot(*main_metric, hotspot_time, DynamicalSystem::LAGRANGIAN, x.at(offset[mode]), x.at(offset[mode] + 1), x.at(offset[mode] + 2) * unit.s, x.at(offset[mode] + 3), false);
	double t = 0., tStep = 0., tRec = 0.25 * unit.s;
	size_t tStepNumber = t1 * 4;
	ublas::bounded_vector<double, 8> position;
	ublas::bounded_vector<double, 5> view_info;
	if (mode == 0) { // circular
		if (hotspot.InitializeCircular(r, x.at(8) * M_PI / 180., x.at(9)) != StatusCode::SUCCESS) {
			PrintlnError("Initialize Circular Error!");
			return HSExit(x);
		}
	} else if (mode == 1) { // geodesic
		if (hotspot.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., x.at(11), x.at(12), inclination, rotation) != StatusCode::SUCCESS) {
			PrintlnError("Initialize Geodesic Error!");
			return HSExit(x);
		}
	} else if (mode == 2) { // conical
		if (hotspot.InitializeConical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11), false) != StatusCode::SUCCESS) {
			PrintlnError("Initialize Conical Error!");
			return HSExit(x);
		}
	} else if (mode == 3) { // cylindrical
		if (hotspot.InitializeCylindrical(r, x.at(8) * M_PI / 180., x.at(9) * R * unit.mas, x.at(10), x.at(11), false) != StatusCode::SUCCESS) {
			PrintlnError("Initialize Cylindrical Error!");
			return HSExit(x);
		}
	} else if (mode == 4) { // conical fixed energy
		if (hotspot.InitializeConical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11), true) != StatusCode::SUCCESS) {
			PrintlnError("Initialize Conical Fixed Energy Error!");
			return HSExit(x);
		}
	} else if (mode == 5) { // cylindrical fixed energy
		if (hotspot.InitializeCylindrical(r, x.at(8) * M_PI / 180., x.at(9) * R * unit.mas, x.at(10), x.at(11), true) != StatusCode::SUCCESS) {
			PrintlnError("Initialize Cylindrical Fixed Energy Error!");
			return HSExit(x);
		}
	} else if (mode == 6) { // schwarzschild geodesic with v_r/v_phi replaced with energy and angular momentum
		PRECISION dtau_dt = (1. - 2. / r) / x.at(11);
		PRECISION dphi_dt = x.at(12) / power2(r) * dtau_dt;
		PRECISION dr_dt = -squareRoot((1. - 2. / r - power2(dtau_dt) - power2(r * dphi_dt)) * (1. - 2. / r));
		if (hotspot.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., dr_dt, dphi_dt, inclination, rotation) != StatusCode::SUCCESS)
			return HSExit(x);
	} else if (mode == 7) { // conical with last dim as angular momentum
		PRECISION theta = x.at(8) * M_PI / 180.;
		if (hotspot.InitializeConical(r, theta, x.at(9) * M_PI / 180., x.at(10), x.at(11) / power2(r * sin(theta)), false) != StatusCode::SUCCESS)
			return HSExit(x);
	} else {
		PrintlnWarning("Mode Error!");
		return py::array_t<double>();
	}
	hotspot.Position(position);
	double h = 1.;
	int status = 0;
	auto result = py::array_t<double>(tStepNumber * HOTSPOT_RETURN_WIDTH);
	double *result_ptr = result.mutable_data();
	ublas::bounded_vector<double, 8> result_position;
	for (size_t i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		if (status = hotspot.IntegratorApply(t, tStep, h); status != 0)
			py::print("[!] IntegratorApply status =", status);
		hotspot.Position(position);
		if (gr_time_delay) {
			if (view_ptr->Trace(position, hotspot_time, view_info, false) != StatusCode::SUCCESS) {
				PrintlnError("Trace star this position Error!");
				result_ptr[0] = NAN;
				result_ptr += HOTSPOT_RETURN_WIDTH;
				continue;
			}
			view_info[3] /= unit.s;
			copy(view_info.begin(), view_info.end(), result_ptr + 10);
			// result_ptr[15] = Flux(hotspot.Luminosity(tStep), result_ptr[14], result_ptr[12]);
			// result_ptr[16] = FluxDensity(hotspot.SpectralDensity(tStep, result_ptr[12]), result_ptr[14]);
		}
		sphericalToCartesian(position);
		copy(position.begin(), position.end(), result_ptr);
		result_ptr[8] = t / unit.s;
		const double delta_epsilon = 1. - 2. / hypot(result_ptr[1], result_ptr[2], result_ptr[3]);
		result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - hypotSquare(result_ptr[5], result_ptr[6], result_ptr[7])) - 1.) * 299792.458;
		result_ptr[15] = hotspot.Energy();
		result_ptr[16] = hotspot.AngularMomentum();
		result_ptr += HOTSPOT_RETURN_WIDTH;
	}
	return result.reshape({tStepNumber, HOTSPOT_RETURN_WIDTH});
}

tuple<py::array_t<double>, int> CalculateHSOrbit(const py::array_t<double> &x, int metric, int mode, bool ray_tracing, bool gr_time_delay, const py::array_t<double> &obs_time, bool calculate_magnification) {
	assert(0 <= mode && mode < 8);
#ifdef SBODY_RELEASE
	auto t_start = chrono::steady_clock::now();
#endif
	const array<double, 9> estimate_step = {120., 60., 30., 10., 3., 1., 0.5, 0.1, 0.};
	size_t estimate_idx = 0;
	PRECISION this_gr_obs_time, last_gr_obs_time, gr_offset = 0.;
	Unit unit(x.at(0)); // double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination,
	PRECISION R = x.at(2) * unit.pc, r = x.at(7) * R * unit.mas, inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
	const array<int, 8> offset = {10, 13, 12, 12, 12, 12, 13, 12};
	unique_ptr<Metric<PRECISION>> main_metric;
	unique_ptr<View<PRECISION>> view_ptr;
	if (metric == 1 || (metric == 2 && abs(x.at(1)) < boost::math::tools::forth_root_epsilon<double>()))
		main_metric = make_unique<Schwarzschild<PRECISION>>();
	else if (metric == 2)
		main_metric = make_unique<Kerr<PRECISION>>(x.at(1));
	else {
		PrintlnError("Metric Error!");
		return {py::array_t<double>(), StatusCode::INVALID_ARGUMENT};
	}
	view_ptr = make_unique<View<PRECISION>>(*main_metric, R, inclination, rotation);
	TimeSystem hotspot_time = TimeSystem::T;
	Hotspot<PRECISION> hotspot(*main_metric, hotspot_time, DynamicalSystem::LAGRANGIAN, x.at(offset[mode]), x.at(offset[mode] + 1), x.at(offset[mode] + 2) * unit.s, x.at(offset[mode] + 3) * unit.s, false);
	PRECISION t = 0., tStep = 0., tRec = 0.1 * unit.s, t0, last_obs_time = 0., this_obs_time;
	ublas::bounded_vector<PRECISION, 8> last_position, this_position, cartesian_position;
	ublas::bounded_vector<PRECISION, 5> last_view_info, this_view_info;
	ublas::bounded_vector<PRECISION, 2> last_hotspot_info, this_hotspot_info;
	const PRECISION sin_inc = sin(inclination), cos_inc = cos(inclination);
	if (mode == 0) { // circular
		if (int status = hotspot.InitializeCircular(r, x.at(8) * M_PI / 180., x.at(9)); status != StatusCode::SUCCESS)
			return {HSExit(x), status};
	} else if (mode == 1) { // geodesic
		if (int status = hotspot.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., x.at(11), x.at(12), inclination, rotation); status != StatusCode::SUCCESS)
			return {HSExit(x), status};
	} else if (mode == 2) { // conical
		if (int status = hotspot.InitializeConical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11), false); status != StatusCode::SUCCESS)
			return {HSExit(x), status};
	} else if (mode == 3) { // cylindrical
		if (int status = hotspot.InitializeCylindrical(r, x.at(8) * M_PI / 180., x.at(9) * R * unit.mas, x.at(10), x.at(11), false); status != StatusCode::SUCCESS)
			return {HSExit(x), status};
	} else if (mode == 4) { // conical fixed energy
		if (int status = hotspot.InitializeConical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11), true); status != StatusCode::SUCCESS)
			return {HSExit(x), status};
	} else if (mode == 5) { // cylindrical fixed energy
		if (int status = hotspot.InitializeCylindrical(r, x.at(8) * M_PI / 180., x.at(9) * R * unit.mas, x.at(10), x.at(11), true); status != StatusCode::SUCCESS)
			return {HSExit(x), status};
	} else if (mode == 6) { // schwarzschild geodesic with v_r/v_phi replaced with energy and angular momentum
		PRECISION dtau_dt = (1. - 2. / r) / x.at(11);
		PRECISION v_phi = x.at(12) / r * dtau_dt;
		PRECISION v_r = -squareRoot((1. - 2. / r - power2(dtau_dt) - power2(v_phi)) * (1. - 2. / r));
		if (int status = hotspot.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., v_r, v_phi, inclination, rotation); status != StatusCode::SUCCESS)
			return {HSExit(x), status};
	} else if (mode == 7) { // conical with last dim as angular momentum
		PRECISION theta = x.at(8) * M_PI / 180.;
		if (int status = hotspot.InitializeConical(r, theta, x.at(9) * M_PI / 180., x.at(10), x.at(11) / power2(r * sin(theta)), false); status != StatusCode::SUCCESS)
			return {HSExit(x), status};
	} else {
		PrintlnWarning("Mode Error!");
		return {py::array_t<double>(), StatusCode::INVALID_ARGUMENT};
	}
	hotspot.Position(last_position);
	if (last_position[1] < 3.) {
		// py::print("[!] R < 3.0, stop!");
		return {py::array_t<double>(), StatusCode::FALL_INTO_BLACK_HOLE};
	}
	PRECISION x0 = last_position[1] * abs(sin(last_position[2])) * cos(last_position[3]), z0 = last_position[1] * cosTheta(last_position[2]);
	if (gr_time_delay) {
		if (int status = view_ptr->Trace(last_position, hotspot_time, last_view_info, calculate_magnification); status != StatusCode::SUCCESS) {
			// PrintlnError("Initially trace star Error!");
			auto numpy = File(FileFormat::NUMPY, fmt::format("log/hotspot/metric_{}/mode_{}/trace/{}", metric, mode, chrono::duration_cast<chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()), {6});
			if (metric == 2)
				numpy.save(x.at(1));
			else
				numpy.save(0.0);
			numpy.save(inclination);
			for (int i = 0; i < 4; ++i)
				numpy.save(last_position[i]);
			return {HSExit(x), status};
		}
		t0 = last_view_info[3] / unit.s;
	}
	PRECISION vz0 = x0 * sin_inc + z0 * cos_inc;
	PRECISION h = -1.;
	int status = 0;
	h = 1.;
	size_t idx = 0, size = obs_time.size();
	auto result = py::array_t<double>(size * HOTSPOT_RETURN_WIDTH);
	double *result_ptr = result.mutable_data();
	for (int i = 0;; ++i) {
#ifdef SBODY_RELEASE
		if (auto now = chrono::steady_clock::now(); now - t_start > chrono::seconds(5)) {
			// py::print("CalculateHSOrbit() timeout!");
			// auto numpy = NumPy(fmt::format("log/hotspot/metric_{}/mode_{}/timeout/{}", metric, mode, chrono::duration_cast<chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()), {static_cast<size_t>(x.size())});
			// for (int j = 0; j < x.size(); ++j)
			// 	numpy.Save(x.at(j));
			return {HSExit(x), StatusCode::COMPUTATION_TIMEOUT};
		}
#endif
		tStep += tRec;
		if (status = hotspot.IntegratorApply(t, tStep, h); status != 0) {
			// py::print("[!] IntegratorApply status =", status);
			return {py::array_t<double>(), status};
		}
#ifndef SBODY_RELEASE
		hotspot.Normalize();
#endif
		hotspot.Position(this_position);
		if (this_position[1] < 3.) {
			// py::print("[!] R < 3.0, stop!");
			return {py::array_t<double>(), StatusCode::FALL_INTO_BLACK_HOLE};
		}
		this_obs_time = (t + vz0 - this_position[1] * (sinTheta(this_position[2]) * cos(this_position[3]) * sin_inc + cosTheta(this_position[2]) * cos_inc)) / unit.s;
		if (gr_time_delay) {
			if (this_obs_time + gr_offset + estimate_step[estimate_idx] > obs_time.at(idx)) {
				if (status = view_ptr->Trace(this_position, hotspot_time, this_view_info, calculate_magnification); status != StatusCode::SUCCESS) {
					if (this_obs_time + gr_offset - 0.1 > obs_time.at(idx)) {
						// PrintlnError("Trace star this position Error!");
						auto numpy = File(FileFormat::NUMPY, fmt::format("log/hotspot/metric_{}/mode_{}/trace/{}", metric, mode, chrono::duration_cast<chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()), {6});
						if (metric == 2)
							numpy.save(x.at(1));
						else
							numpy.save(0.0);
						numpy.save(inclination);
						for (int j = 0; j < 4; ++j)
							numpy.save(this_position[j]);
						return {HSExit(x), status};
					}
					if (estimate_idx + 2 < estimate_step.size())
						++estimate_idx;
					continue;
				}
				this_view_info[3] /= unit.s;
				this_gr_obs_time = t / unit.s + t0 - this_view_info[3];
				// if (calculate_magnification) {
				// 	this_view_info[5] = Flux(hotspot.Luminosity(tStep), this_view_info[4], this_view_info[2]);
				// 	this_view_info[6] = FluxDensity(hotspot.SpectralDensity(tStep, this_view_info[2]), this_view_info[4]);
				// }
				this_hotspot_info[0] = hotspot.Energy();
				this_hotspot_info[1] = hotspot.AngularMomentum();
				gr_offset = this_gr_obs_time - this_obs_time;
				if (abs(gr_offset) > 3600.) { // difference between Romer delay and full GR delay
					// PrintlnError("gr_offset > 3600 s!");
					// for (int j = 0; j < 8; ++j)
					// 	py::print("this_position[", j, "]=", this_position[j]);
					auto numpy = File(FileFormat::NUMPY, fmt::format("log/hotspot/metric_{}/mode_{}/gr_offset/{}", metric, mode, chrono::duration_cast<chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()), {6});
					if (metric == 2)
						numpy.save(x.at(1));
					else
						numpy.save(0.0);
					numpy.save(inclination);
					for (int j = 0; j < 4; ++j)
						numpy.save(this_position[j]);
					return {HSExit(x), StatusCode::FAILURE};
				}
				if (++estimate_idx >= estimate_step.size()) {
					estimate_idx = 0;
					if (min(abs(obs_time.at(idx) - last_gr_obs_time), abs(obs_time.at(idx) - this_gr_obs_time)) > 10. * abs(this_gr_obs_time - last_gr_obs_time))
						return {HSExit(x), StatusCode::NUMERIC_ERROR};
					InterpolateSphericalPositionToCartesian(static_cast<PRECISION>(obs_time.at(idx)), last_gr_obs_time, this_gr_obs_time, last_position, this_position, cartesian_position);
					copy(cartesian_position.begin(), cartesian_position.end(), result_ptr);
					result_ptr[8] = LinearInterpolation(static_cast<PRECISION>(obs_time.at(idx)), last_gr_obs_time, this_gr_obs_time, tStep - tRec, tStep) / unit.s;
					const PRECISION delta_epsilon = 1. - 2. / hypot(result_ptr[1], result_ptr[2], result_ptr[3]);
					result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - hypotSquare(result_ptr[5], result_ptr[6], result_ptr[7])) - 1.) * 299792.458;
					LinearInterpolation(static_cast<PRECISION>(obs_time.at(idx)), last_gr_obs_time, this_gr_obs_time, last_view_info, this_view_info, result_ptr + 10);
					LinearInterpolation(static_cast<PRECISION>(obs_time.at(idx)), last_gr_obs_time, this_gr_obs_time, last_hotspot_info, this_hotspot_info, result_ptr + 15);
					result_ptr += HOTSPOT_RETURN_WIDTH;
					if (++idx >= size)
						break;
				}
				last_gr_obs_time = this_gr_obs_time;
				swap(last_view_info, this_view_info);
				swap(last_hotspot_info, this_hotspot_info);
			}
		} else if (this_obs_time > obs_time.at(idx)) {
			InterpolateSphericalPositionToCartesian(static_cast<PRECISION>(obs_time.at(idx)), last_obs_time, this_obs_time, last_position, this_position, cartesian_position);
			copy(cartesian_position.begin(), cartesian_position.end(), result_ptr);
			if (ray_tracing) {
				cartesianToSpherical(cartesian_position);
				if (status = view_ptr->Trace(cartesian_position, hotspot_time, last_view_info, calculate_magnification); status != StatusCode::SUCCESS) {
					PrintlnError("Trace star position Error!");
					return {HSExit(x), status};
				}
				last_view_info[3] /= unit.s;
				copy(last_view_info.begin(), last_view_info.end(), result_ptr + 10);
				// if (calculate_magnification) {
				// 	result_ptr[15] = Flux(hotspot.Luminosity(tStep), static_cast<PRECISION>(result_ptr[14]), static_cast<PRECISION>(result_ptr[12]));
				// 	result_ptr[16] = FluxDensity(hotspot.SpectralDensity(tStep, result_ptr[12]), static_cast<PRECISION>(result_ptr[14]));
				// }
				result_ptr[15] = hotspot.Energy();
				result_ptr[16] = hotspot.AngularMomentum();
			}
			result_ptr[8] = LinearInterpolation(static_cast<PRECISION>(obs_time.at(idx)), last_obs_time, this_obs_time, tStep - tRec, tStep) / unit.s;
			const double delta_epsilon = 1. - 2. / hypot(result_ptr[1], result_ptr[2], result_ptr[3]);
			result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) / sqrt(delta_epsilon - hypotSquare(result_ptr[5], result_ptr[6], result_ptr[7])) - 1.) * 299792.458;
			result_ptr += HOTSPOT_RETURN_WIDTH;
			if (++idx >= size)
				break;
		}
		last_obs_time = this_obs_time;
		swap(last_position, this_position);
	}
	return {result.reshape({size, HOTSPOT_RETURN_WIDTH}), 0};
}

py::array_t<double> CalculateHotspotEnergyAngularMomentumAndMinimumRadius(const py::array_t<double> &x, int metric) {
	Unit unit(x.at(0)); // double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination,
	double R = x.at(2) * unit.pc, r = x.at(7) * R * unit.mas, inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
	unique_ptr<Metric<double>> main_metric;
	auto result = py::array_t<double>(3);
	double *result_ptr = result.mutable_data();
	result_ptr[0] = result_ptr[1] = NAN;
	result_ptr[2] = 0.0;
	if (metric == 1 || (metric == 2 && abs(x.at(1)) < boost::math::tools::forth_root_epsilon<double>()))
		main_metric = make_unique<Schwarzschild<double>>();
	else if (metric == 2)
		main_metric = make_unique<Kerr<double>>(x.at(1));
	else {
		PrintlnError("Metric Error!");
		return result;
	}
	Hotspot<double> hotspot(*main_metric, TimeSystem::T, DynamicalSystem::LAGRANGIAN, x.at(13), x.at(14), x.at(15) * unit.s, x.at(16) * unit.s, false);
	if (hotspot.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., x.at(11), x.at(12), inclination, rotation) != StatusCode::SUCCESS)
		return result;
	result_ptr[0] = hotspot.Energy();
	if (metric == 1) {
		result_ptr[1] = power2(hotspot.Position(1)) * sqrt(power2(hotspot.Position(6)) + power2(std::sin(hotspot.Position(2)) * hotspot.Position(7)));
		double cubic_roots[3];
		double coeff = 1. / (power2(result_ptr[0]) - 1.);
		double L2_coeff = power2(result_ptr[1]) * coeff;
		int num_of_roots = polySolveCubic(2. * coeff, -L2_coeff, 2. * L2_coeff, cubic_roots);
		for (int i = 0; i < num_of_roots; ++i)
			if (cubic_roots[i] > 0.) {
				result_ptr[2] = cubic_roots[i];
				break;
			}
	} else
		result_ptr[1] = hotspot.AngularMomentum();
	return result;
}

double HSFluxChi2(py::array_t<double> &x, py::array_t<double> &intrinsic_time, py::array_t<double> &obs_magnification, py::array_t<double> &obs_redshift, py::array_t<double> &obs_flux, py::array_t<double> &flux_sigma) {
	Unit unit(x.at(0));
	const int size = obs_flux.size();
	unique_ptr<Metric<double>> main_metric;
	Hotspot<double> hotspot(*main_metric, TimeSystem::TAU, DynamicalSystem::LAGRANGIAN, x.at(0), x.at(1), x.at(2) * unit.s, x.at(3) * unit.s, false);
	double chi2 = 0.;
	for (int i = 0; i < size; ++i)
		chi2 += power2((FluxDensity(hotspot.SpectralDensity(intrinsic_time.at(i), obs_redshift.at(i)), obs_magnification.at(i)) - obs_flux.at(i)) / flux_sigma.at(i));
	return chi2;
}
tuple<py::array_t<double>, int> HSDistinctChi2(py::array_t<double> &x, int metric, int mode, int gr_switch, py::array_t<double> &obs_time, py::array_t<double> &obs_flux, py::array_t<double> &flux_sigma, py::array_t<double> &obs_ra, py::array_t<double> &ra_sigma, py::array_t<double> &obs_dec, py::array_t<double> &dec_sigma, double disk_flux) {
	if (metric == 0 && gr_switch > 0) {
		x.mutable_at(1) = 1.;
		gr_switch = 0;
	}
	auto [obs_data, status] = CalculateHSOrbit(x, metric, mode, gr_switch > 0, (gr_switch & 4) > 0, obs_time, gr_switch & 8);
	if (obs_data.size() == 0)
		return {py::array_t<double>(), status};
	const int size = obs_flux.size();
	auto result = py::array_t<double>(size);
	double *result_ptr = result.mutable_data();
	memset(result_ptr, 0, size * sizeof(double));
	if (gr_switch & 8)
		for (int i = 0; i < size; ++i)
			result_ptr[i] += power2((obs_flux.at(i) - obs_data.at(i, 16)) / flux_sigma.at(i));
	if (gr_switch & 2)
		for (int i = 0; i < size; ++i) { // 9.870628713769018e-6 = 1. / (unit.pc * unit.mas)
			result_ptr[i] += power2((obs_ra.at(i) + obs_data.at(i, 10) * x.at(0) * 9.870628713769018e-6 / x.at(2) * (1. - disk_flux / obs_flux.at(i)) - x.at(3)) / ra_sigma.at(i));
			result_ptr[i] += power2((obs_dec.at(i) - obs_data.at(i, 11) * x.at(0) * 9.870628713769018e-6 / x.at(2) * (1. - disk_flux / obs_flux.at(i)) - x.at(4)) / dec_sigma.at(i));
		}
	else {
		PrintlnError("Inclination and rotation is ignored!");
		for (int i = 0; i < size; ++i) {
			result_ptr[i] += power2((obs_ra.at(i) + obs_data.at(i, 2) * x.at(0) * 9.870628713769018e-6 / x.at(2) * (1. - disk_flux / obs_flux.at(i)) - x.at(3)) / ra_sigma.at(i));	 // FIXME: inclination!
			result_ptr[i] += power2((obs_dec.at(i) + obs_data.at(i, 1) * x.at(0) * 9.870628713769018e-6 / x.at(2) * (1. - disk_flux / obs_flux.at(i)) - x.at(4)) / dec_sigma.at(i)); // FIXME: inclination!
		}
	}
	return {result, 0};
}
py::tuple HotspotChi2(py::array_t<double> &x, int metric, int mode, int gr_switch, py::array_t<double> &obs_time, py::array_t<double> &obs_flux, py::array_t<double> &flux_sigma, py::array_t<double> &obs_ra, py::array_t<double> &ra_sigma, py::array_t<double> &obs_dec, py::array_t<double> &dec_sigma, double black_hole_sigma) {
	cerr.rdbuf(cout.rdbuf());
	auto [result, status] = HSDistinctChi2(x, metric, mode, gr_switch, obs_time, obs_flux, flux_sigma, obs_ra, ra_sigma, obs_dec, dec_sigma, 0.0);
	if (result.size() == 0)
		return py::make_tuple(std::numeric_limits<double>::infinity(), status);
	double prob_sum = 0.0;
	for (int i = 0; i < result.size(); ++i)
		prob_sum += result.at(i);
	if (gr_switch & 16)
		prob_sum += (power2(x.at(3)) + power2(x.at(4))) / power2(black_hole_sigma);
	if (isnan(prob_sum)) {
		for (int i = 0; i < result.size(); ++i)
			py::print(result.at(i));
		return py::make_tuple(std::numeric_limits<double>::infinity(), -1);
	}
	return py::make_tuple(prob_sum, 0);
}

py::array_t<double> FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_target, const double theta_target, const double phi_target) {
	ublas::bounded_vector<double, 8> photon;
	double alpha = 0., beta = 0., photon_time;
	Schwarzschild<double>().FastTrace(r_observer, theta_observer, sin_theta_observer, cos_theta_observer, r_target, theta_target, phi_target, alpha, beta, photon, photon_time);
	auto res = py::array_t<double>(9);
	copy(photon.begin(), photon.end(), res.mutable_data());
	res.mutable_data()[8] = photon_time;
	return res;
}

py::array_t<double> CppSort(py::array_t<double> data) {
	double *data_ptr = data.mutable_data();
	sort(data_ptr, data_ptr + data.size());
	return data;
}

py::array_t<double> KerrCoverage(double spin, double theta) {
	namespace ublas = boost::numeric::ublas;
	Unit unit(4.15e6);
	auto main_metric = Kerr<double>(spin);
	const double sin_theta = sin(theta), cos_theta = cos(theta);
	ublas::bounded_vector<double, 201> u;
	ProgressBar bars(2, "Kerr Coverage");
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
			main_metric.CalcThetaPhi(8.18 * unit.pc, theta, sin_theta, cos_theta, alpha, beta, u, result.mutable_data() + offset_0, result.mutable_data() + offset_1, result.mutable_data() + offset_2, result.mutable_data() + offset_3);
			offset_0 += offset;
			offset_1 += offset;
			offset_2 += offset;
			offset_3 += offset;
			bars[0].set_progress((i + 100 + (j + 100) / 201.) * 100. / 201.);
		}
	}
	return result.reshape({201, 201, 4, 201});
}

py::array_t<double> KerrRingShape(const py::array_t<double> &r, double R, double spin, double theta) {
	namespace ublas = boost::numeric::ublas;
	Unit unit(4.15e6);
	auto main_metric = Kerr<double>(spin);
	const double sin_theta = sin(theta), cos_theta = cos(theta);
	int N = r.size();
	ublas::vector<double> u(N);
	auto result = py::array_t<double>(101 * 4 * N);
	memset(result.mutable_data(), 0xff, sizeof(double) * 101 * 4 * N);
	if (!isnan(result.data()[0]))
		return py::array_t<double>();
	for (int i = 0; i < N; ++i)
		u[i] = 1. / r.at(i);
	int offset = 4 * N;
	int offset_0 = 0;
	int offset_1 = N;
	int offset_2 = offset_1 * 2;
	int offset_3 = offset_1 * 3;
	for (int i = 0; i < 101; ++i) {
		double alpha_beta_theta = boost::math::constants::two_pi<double>() * i / 101.;
		double alpha = R * cos(alpha_beta_theta);
		double beta = R * sin(alpha_beta_theta);
		main_metric.CalcThetaPhi(8.18 * unit.pc, theta, sin_theta, cos_theta, alpha, beta, u, result.mutable_data() + offset_0, result.mutable_data() + offset_1, result.mutable_data() + offset_2, result.mutable_data() + offset_3);
		offset_0 += offset;
		offset_1 += offset;
		offset_2 += offset;
		offset_3 += offset;
	}
	return result.reshape({101, 4, N});
}

py::array_t<double> CalculateShadowSilhouette(const py::array_t<double> &black_hole_parameters, int metric) {
	return py::array_t<double>();
}

py::array_t<double> ShadowChi2(const py::array_t<double> &observed_shadow_silhouette, const py::array_t<double> &black_hole_parameters, int metric, double resolution) {
	int data_point_number = observed_shadow_silhouette.shape(0);
	auto observed_shadow_silhouette_reader = observed_shadow_silhouette.unchecked<2>();
	vector<double> corresponding_weight(data_point_number);
	resolution = 0.5 / resolution;
	for (int i = 1; i < data_point_number; ++i) {
		double half_distance = resolution * sqrt(power2(observed_shadow_silhouette_reader(i, 0) - observed_shadow_silhouette_reader(i - 1, 0)) + power2(observed_shadow_silhouette_reader(i, 1) - observed_shadow_silhouette_reader(i - 1, 1)));
		corresponding_weight[i] += half_distance;
		corresponding_weight[i - 1] += half_distance;
	}
	double half_distance = resolution * sqrt(power2(observed_shadow_silhouette_reader(data_point_number - 1, 0) - observed_shadow_silhouette_reader(0, 0)) + power2(observed_shadow_silhouette_reader(data_point_number - 1, 1) - observed_shadow_silhouette_reader(0, 1)));
	corresponding_weight[data_point_number - 1] += half_distance;
	corresponding_weight[0] += half_distance;
	auto result = py::array_t<double>(data_point_number);
	auto result_writer = result.mutable_unchecked();
	auto calculated_shadow_silhouette = CalculateShadowSilhouette(black_hole_parameters, metric);
	auto calculated_shadow_silhouette_reader = calculated_shadow_silhouette.unchecked<2>();
	for (int i = 0; i < data_point_number; ++i) {
		result_writer(i) = corresponding_weight[i] * sqrt(power2(calculated_shadow_silhouette_reader(i, 0) - observed_shadow_silhouette_reader(i, 0)) + power2(calculated_shadow_silhouette_reader(i, 1) - observed_shadow_silhouette_reader(i, 1)));
	}
	return result;
}

int HotspotOrbitTest(const py::array_t<double> &x, int metric, int mode, bool ray_tracing, bool gr_time_delay) {
	assert(0 <= mode && mode < 6);
	Unit unit(x.at(0)); // double fSP, double R, double r, double theta, double phi, double v_r, double v_phi, double inclination,
	PRECISION R = x.at(2) * unit.pc, r = x.at(7) * R * unit.mas, inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
	const array<int, 6> offset = {10, 13, 12, 12, 12, 12};
	unique_ptr<Metric<PRECISION>> main_metric;
	unique_ptr<View<PRECISION>> view_ptr;
	if (metric == 1 || (metric == 2 && abs(x.at(1)) < boost::math::tools::forth_root_epsilon<double>()))
		main_metric = make_unique<Schwarzschild<PRECISION>>();
	else if (metric == 2)
		main_metric = make_unique<Kerr<PRECISION>>(x.at(1));
	else {
		PrintlnError("Metric Error!");
		return StatusCode::FAILURE;
	}
	int status = StatusCode::SUCCESS;
	TimeSystem hotspot_time = TimeSystem::T;
	Hotspot<PRECISION> hotspot(*main_metric, hotspot_time, DynamicalSystem::LAGRANGIAN, x.at(offset[mode]), x.at(offset[mode] + 1), x.at(offset[mode] + 2) * unit.s, x.at(offset[mode] + 3) * unit.s, false);
	PRECISION t = 0., tStep = 0., tRec = 10. * unit.s;
	if (mode == 0) { // circular
		if (hotspot.InitializeCircular(r, x.at(8) * M_PI / 180., x.at(9)) != StatusCode::SUCCESS)
			return StatusCode::FAILURE;
	} else if (mode == 1) { // geodesic
		if (hotspot.InitializeGeodesic(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10) * M_PI / 180., x.at(11), x.at(12), inclination, rotation) != StatusCode::SUCCESS)
			return StatusCode::FAILURE;
	} else if (mode == 2) { // conical
		if (hotspot.InitializeConical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11), false) != StatusCode::SUCCESS)
			return StatusCode::FAILURE;
	} else if (mode == 3) { // cylindrical
		if (hotspot.InitializeCylindrical(r, x.at(8) * M_PI / 180., x.at(9) * R * unit.mas, x.at(10), x.at(11), false) != StatusCode::SUCCESS)
			return StatusCode::FAILURE;
	} else if (mode == 4) { // conical fixed energy
		if (hotspot.InitializeConical(r, x.at(8) * M_PI / 180., x.at(9) * M_PI / 180., x.at(10), x.at(11), true) != StatusCode::SUCCESS)
			return StatusCode::FAILURE;
	} else if (mode == 5) { // cylindrical fixed energy
		if (hotspot.InitializeCylindrical(r, x.at(8) * M_PI / 180., x.at(9) * R * unit.mas, x.at(10), x.at(11), true) != StatusCode::SUCCESS)
			return StatusCode::FAILURE;
	} else {
		PrintlnWarning("Mode Error!");
		return StatusCode::FAILURE;
	}
	PRECISION h = 1.0;
	for (int i = 0; i < 100; ++i) {
		tStep += tRec;
		if (status = hotspot.IntegratorApply(t, tStep, h); status != 0) {
			py::print("[!] IntegratorApply status =", status);
			return StatusCode::FAILURE;
		}
		if (hotspot.Position(1) < 3.)
			break;
#ifndef SBODY_RELEASE
		hotspot.Normalize();
#endif
		py::print("r:", hotspot.Position(1), "Energy:", hotspot.Energy(), "Angular Momentum:", hotspot.AngularMomentum());
	}
	return StatusCode::SUCCESS;
}

void bindHotspotModule(py::module_ &m) {
	m.def("CalculateFullHSOrbit", &CalculateFullHSOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("CalculateHSOrbit", &CalculateHSOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

	m.def("HotspotChi2", &HotspotChi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

	m.def("HSDistinctChi2", &HSDistinctChi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

	m.def("CalculateHotspotEnergyAngularMomentumAndMinimumRadius", &CalculateHotspotEnergyAngularMomentumAndMinimumRadius, R"pbdoc(
		Calculate Hotspot Energy Angular Momentum And Minimum Radius for geodesic orbits.
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

	m.def("KerrRingShape", &KerrRingShape, R"pbdoc(
		KerrRingShape
	)pbdoc");

	m.def("HotspotOrbitTest", &HotspotOrbitTest, R"pbdoc(
		HotspotOrbitTest
	)pbdoc");
}
