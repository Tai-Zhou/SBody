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
#include <filesystem>
#include <memory>
#include <pybind11/cast.h>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <boost/algorithm/algorithm.hpp>
#include <fmt/core.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "file/file.hpp"
#include "infra/infra_context.hpp"
#include "infra/status.hpp"
#include "math/geometry.hpp"
#include "math/interpolation.hpp"
#include "math/precision.hpp"
#include "math/solver/polynomial.hpp"
#include "metric/metrics.hpp"
#include "object/particle.hpp"
#include "units.hpp"
#include "view/flux.hpp"
#include "view/view.hpp"

namespace py = pybind11;
namespace ublas = boost::numeric::ublas;
using namespace sbody;

constexpr real_t HOTSPOT_TRACE_INNER_RADIUS = 3.0;
constexpr real_t M_SUN_RG_IN_PARSEC_MAS = 1.0 / (20896825544594.496 * mas_in_rad);
constexpr real_t SPEED_OF_LIGHT_IN_KM_PER_SECOND = 299792.458;

constexpr size_t HOTSPOT_RETURN_WIDTH = 17UL;

constexpr size_t OBS_TIME_IDX = 0;
constexpr size_t OBS_FLUX_IDX = 1;
constexpr size_t OBS_FLUX_SIGMA_IDX = 2;
constexpr size_t OBS_RA_IDX = 3;
constexpr size_t OBS_RA_SIGMA_IDX = 4;
constexpr size_t OBS_DEC_IDX = 5;
constexpr size_t OBS_DEC_SIGMA_IDX = 6;

constexpr real_t HOTSPOT_TRACE_ABS_ACC = 1e-12;
constexpr real_t HOTSPOT_TRACE_REL_ACC = 1e-12;
constexpr int TIME_OUT_SECONDS = 5;
constexpr LogLevel HOTSPOT_LOG_LEVEL = LogLevel::ERROR;
constexpr FILE *HOTSPOT_STD_OUT = nullptr;
constexpr FILE *HOTSPOT_STD_ERR = nullptr;

std::vector<real_t> convertFromNumpyDoubleArray(py::array_t<double> array);

std::vector<uint8_t> convertFromNumpyBoolArray(py::array_t<bool> array);

py::array_t<double> hotspotExit(const std::vector<real_t> &x, const InfraContext &infra_context) {
    for (size_t i = 0; i < x.size(); ++i)
        infra_context.logger_.info("x[{}]={}", i, x[i]);
    return py::array_t<double>();
}

std::unique_ptr<Metric<real_t>> getMetric(int metric, real_t spin) {
    if (metric == 1 || (metric == 2 && abs(spin) < boost::math::tools::forth_root_epsilon<real_t>()))
        return std::make_unique<Schwarzschild<real_t>>();
    else if (metric == 2)
        return std::make_unique<Kerr<real_t>>(spin);
    throw std::invalid_argument("Metric Error!");
}

int initializeHotspotMotion(const Unit<real_t> &unit, Hotspot<real_t> &hotspot, const std::vector<real_t> &x, int mode) {
    const real_t R = x[2] * unit.pc;
    const real_t r = x[7] * R * mas_in_rad;
    const real_t inclination = x[5] * M_PI / 180.;
    const real_t rotation = x[6] * M_PI / 180.;
    if (mode == 0) { // circular
        return hotspot.initializeCircular(r, x[8] * M_PI / 180., x[9]);
    } else if (mode == 1) { // geodesic
        return hotspot.initializeGeodesic(r, x[8] * M_PI / 180., x[9] * M_PI / 180., x[10] * M_PI / 180., x[11], x[12],
                                          inclination, rotation);
    } else if (mode == 2) { // conical
        return hotspot.initializeConical(r, x[8] * M_PI / 180., x[9] * M_PI / 180., x[10], x[11], false);
    } else if (mode == 3) { // cylindrical
        return hotspot.initializeCylindrical(r, x[8] * M_PI / 180., x[9] * R * mas_in_rad, x[10], x[11], false);
    } else if (mode == 4) { // conical fixed energy
        return hotspot.initializeConical(r, x[8] * M_PI / 180., x[9] * M_PI / 180., x[10], x[11], true);
    } else if (mode == 5) { // cylindrical fixed energy
        return hotspot.initializeCylindrical(r, x[8] * M_PI / 180., x[9] * R * mas_in_rad, x[10], x[11], true);
    } else if (mode == 6) { // schwarzschild geodesic with v_r/v_phi replaced with energy and angular momentum
        real_t dtau_dt = (1. - 2. / r) / x[11];
        real_t dphi_dt = x[12] / power2(r) * dtau_dt;
        real_t dr_dt = -squareRoot((1. - 2. / r - power2(dtau_dt) - power2(r * dphi_dt)) * (1. - 2. / r));
        return hotspot.initializeGeodesic(r, x[8] * M_PI / 180., x[9] * M_PI / 180., x[10] * M_PI / 180., dr_dt, dphi_dt,
                                          inclination, rotation);
    } else if (mode == 7) { // conical with last dim as angular momentum
        real_t theta = x[8] * M_PI / 180.;
        return hotspot.initializeConical(r, theta, x[9] * M_PI / 180., x[10], x[11] / power2(r * sin(theta)), false);
    }
    return StatusCode::INVALID_ARGUMENT;
}

Hotspot<real_t> initializeHotspot(const Unit<real_t> &unit, Metric<real_t> &metric, const std::vector<real_t> &x,
                                  int motion_mode, const InfraContext &infra_context) {
    const std::array<int, 8> offset = {10, 13, 12, 12, 12, 12, 13, 12};
    Hotspot<real_t> hotspot(metric, TimeSystem::T, DynamicalSystem::LAGRANGIAN, HOTSPOT_TRACE_ABS_ACC, HOTSPOT_TRACE_REL_ACC,
                            x[offset[motion_mode]], x[offset[motion_mode] + 1], x[offset[motion_mode] + 2] * unit.s,
                            x[offset[motion_mode] + 3] * unit.s, false);
    if (motion_mode >= 0 && initializeHotspotMotion(unit, hotspot, x, motion_mode) != StatusCode::SUCCESS) {
        infra_context.logger_.error("initializeHotspot error!");
        throw std::runtime_error("initializeHotspot error!");
    }
    return hotspot;
}

py::array_t<real_t> calculateHotspotFullOrbit(const std::vector<real_t> &x, int metric, int mode, bool gr_time_delay,
                                              real_t t1) {
    Unit unit(x[0]);
    InfraContext infra_context(HOTSPOT_LOG_LEVEL, HOTSPOT_STD_OUT, HOTSPOT_STD_ERR);
    const real_t R = x[2] * unit.pc, inclination = x[5] * M_PI / 180., rotation = x[6] * M_PI / 180.;
    std::unique_ptr<Metric<real_t>> main_metric;
    try {
        main_metric = getMetric(metric, x[1]);
    } catch (const std::invalid_argument &e) {
        infra_context.logger_.error(e.what());
        return py::array_t<double>();
    }
    std::unique_ptr<View<real_t>> view_ptr = std::make_unique<View<real_t>>(*main_metric, R, inclination, rotation);
    Hotspot hotspot = initializeHotspot(unit, *main_metric, x, mode, infra_context);
    // (*main_metric, TimeSystem::T, DynamicalSystem::LAGRANGIAN, HOTSPOT_TRACE_ABS_ACC, HOTSPOT_TRACE_REL_ACC,
    // x.at(offset[mode]), x.at(offset[mode] + 1), x.at(offset[mode] + 2) * unit.s, x.at(offset[mode] + 3), false);
    real_t t = 0., tStep = 0., tRec = 0.25 * unit.s;
    size_t tStepNumber = t1 * 4;
    ublas::bounded_vector<real_t, 8> position;
    ublas::bounded_vector<real_t, 5> view_info;
    hotspot.position(position);
    real_t h = 1.;
    int status = 0;
    auto result = py::array_t<double>(tStepNumber * HOTSPOT_RETURN_WIDTH);
    double *result_ptr = result.mutable_data();
    for (size_t i = 0; i < tStepNumber; ++i) {
        tStep += tRec;
        if (status = hotspot.integratorApply(t, tStep, h); status != StatusCode::SUCCESS)
            infra_context.logger_.warning("[!] IntegratorApply status =", status);
        hotspot.position(position);
        if (gr_time_delay) {
            if (view_ptr->Trace(position, TimeSystem::T, view_info, HOTSPOT_TRACE_ABS_ACC, HOTSPOT_TRACE_REL_ACC, false, true,
                                infra_context) != StatusCode::SUCCESS) {
                infra_context.logger_.error("Trace star this position Error!");
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

        const real_t delta_epsilon = 1. - 2. / hypot(result_ptr[1], result_ptr[2], result_ptr[3]);
        result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) /
                             sqrt(delta_epsilon - hypotSquare(result_ptr[5], result_ptr[6], result_ptr[7])) -
                         1.) *
                        SPEED_OF_LIGHT_IN_KM_PER_SECOND;
        result_ptr[15] = hotspot.energy();
        result_ptr[16] = hotspot.angularMomentum();
        result_ptr += HOTSPOT_RETURN_WIDTH;
    }
    return result.reshape({tStepNumber, HOTSPOT_RETURN_WIDTH});
}

int calculateHotspotOrbitWithGRDelay(
    const Unit<real_t> &unit, View<real_t> &view, const Hotspot<real_t> &hotspot, const std::array<real_t, 9> &estimate_step,
    size_t &estimate_idx, const std::vector<real_t> &obs_time, size_t &idx, real_t this_obs_time, real_t &gr_offset, real_t t,
    real_t t0, real_t t_step, real_t t_rec, const ublas::bounded_vector<real_t, 8> &last_position,
    const ublas::bounded_vector<real_t, 8> &this_position, ublas::bounded_vector<real_t, 8> &cartesian_position,
    ublas::bounded_vector<real_t, 5> &last_view_info, ublas::bounded_vector<real_t, 5> &this_view_info,
    ublas::bounded_vector<real_t, 2> &last_hotspot_info, ublas::bounded_vector<real_t, 2> &this_hotspot_info,
    real_t &last_gr_obs_time, real_t &this_gr_obs_time, double *&result_ptr, bool calculate_magnification,
    const InfraContext &infra_context) {
    if (int status = view.Trace(this_position, TimeSystem::T, this_view_info, HOTSPOT_TRACE_ABS_ACC, HOTSPOT_TRACE_REL_ACC,
                                calculate_magnification, false, infra_context);
        status != StatusCode::SUCCESS) {
        if (this_obs_time + gr_offset - 0.1 > obs_time.at(idx)) {
            return status;
        }
        if (estimate_idx + 2 < estimate_step.size())
            ++estimate_idx;
        return StatusCode::SUCCESS;
    }
    this_view_info[3] /= unit.s;
    this_gr_obs_time = t / unit.s + t0 - this_view_info[3];
    // if (calculate_magnification) {
    // 	this_view_info[5] = Flux(hotspot.Luminosity(tStep), this_view_info[4], this_view_info[2]);
    // 	this_view_info[6] = FluxDensity(hotspot.SpectralDensity(tStep, this_view_info[2]), this_view_info[4]);
    // }
    this_hotspot_info[0] = hotspot.energy();
    this_hotspot_info[1] = hotspot.angularMomentum();
    gr_offset = this_gr_obs_time - this_obs_time;
    if (abs(gr_offset) > 3600.) { // difference between Romer delay and full GR delay
        return StatusCode::FAILURE;
    }
    if (++estimate_idx >= estimate_step.size()) {
        estimate_idx = 0;
        if (std::min(std::abs(obs_time.at(idx) - last_gr_obs_time), std::abs(obs_time.at(idx) - this_gr_obs_time)) >
            10. * std::abs(this_gr_obs_time - last_gr_obs_time))
            return StatusCode::NUMERIC_ERROR;
        interpolateSphericalPositionToCartesian(static_cast<real_t>(obs_time.at(idx)), last_gr_obs_time, this_gr_obs_time,
                                                last_position, this_position, cartesian_position);
        copy(cartesian_position.begin(), cartesian_position.end(), result_ptr);
        result_ptr[8] = linearInterpolation(static_cast<real_t>(obs_time.at(idx)), last_gr_obs_time, this_gr_obs_time,
                                            t_step - t_rec, t_step) /
                        unit.s;
        const real_t delta_epsilon = 1. - 2. / hypot(result_ptr[1], result_ptr[2], result_ptr[3]);
        result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) /
                             sqrt(delta_epsilon - hypotSquare(result_ptr[5], result_ptr[6], result_ptr[7])) -
                         1.) *
                        SPEED_OF_LIGHT_IN_KM_PER_SECOND;
        linearInterpolation(static_cast<real_t>(obs_time.at(idx)), last_gr_obs_time, this_gr_obs_time, last_view_info,
                            this_view_info, result_ptr + 10);
        linearInterpolation(static_cast<real_t>(obs_time.at(idx)), last_gr_obs_time, this_gr_obs_time, last_hotspot_info,
                            this_hotspot_info, result_ptr + 15);
        result_ptr += HOTSPOT_RETURN_WIDTH;
        ++idx;
    }
    last_gr_obs_time = this_gr_obs_time;
    swap(last_view_info, this_view_info);
    swap(last_hotspot_info, this_hotspot_info);
    return StatusCode::SUCCESS;
}

int calculateHotspotOrbitWithRomerDelay(const Unit<real_t> &unit, View<real_t> &view, const Hotspot<real_t> &hotspot,
                                        const std::vector<real_t> &obs_time, size_t &idx, real_t last_obs_time,
                                        real_t this_obs_time, real_t t_step, real_t t_rec,
                                        ublas::bounded_vector<real_t, 8> &last_position,
                                        ublas::bounded_vector<real_t, 8> &this_position,
                                        ublas::bounded_vector<real_t, 8> &cartesian_position,
                                        ublas::bounded_vector<real_t, 5> &last_view_info, double *&result_ptr, bool ray_tracing,
                                        bool calculate_magnification, const InfraContext &infra_context) {
    interpolateSphericalPositionToCartesian(static_cast<real_t>(obs_time.at(idx)), last_obs_time, this_obs_time, last_position,
                                            this_position, cartesian_position);
    copy(cartesian_position.begin(), cartesian_position.end(), result_ptr);
    if (ray_tracing) {
        cartesianToSpherical(cartesian_position);
        if (int status = view.Trace(cartesian_position, TimeSystem::T, last_view_info, HOTSPOT_TRACE_ABS_ACC,
                                    HOTSPOT_TRACE_REL_ACC, calculate_magnification, false, infra_context);
            status != StatusCode::SUCCESS) {
            infra_context.logger_.error("Trace star position Error!");
            return status;
        }
        last_view_info[3] /= unit.s;
        copy(last_view_info.begin(), last_view_info.end(), result_ptr + 10);
        // if (calculate_magnification) {
        // 	result_ptr[15] = Flux(hotspot.Luminosity(tStep), static_cast<real_t>(result_ptr[14]),
        // static_cast<real_t>(result_ptr[12])); 	result_ptr[16] = FluxDensity(hotspot.SpectralDensity(tStep,
        // result_ptr[12]), static_cast<real_t>(result_ptr[14]));
        // }
        result_ptr[15] = hotspot.energy();
        result_ptr[16] = hotspot.angularMomentum();
    }
    result_ptr[8] =
        linearInterpolation(static_cast<real_t>(obs_time.at(idx)), last_obs_time, this_obs_time, t_step - t_rec, t_step) /
        unit.s;
    const real_t delta_epsilon = 1. - 2. / hypot(result_ptr[1], result_ptr[2], result_ptr[3]);
    result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) /
                         sqrt(delta_epsilon - hypotSquare(result_ptr[5], result_ptr[6], result_ptr[7])) -
                     1.) *
                    SPEED_OF_LIGHT_IN_KM_PER_SECOND;
    result_ptr += HOTSPOT_RETURN_WIDTH;
    ++idx;
    return StatusCode::SUCCESS;
}

std::tuple<py::array_t<double>, int> calculateHotspotOrbit(const std::vector<real_t> &x, int metric, int mode, bool ray_tracing,
                                                           bool gr_time_delay, const std::vector<real_t> &obs_time,
                                                           bool calculate_magnification, const InfraContext &infra_context) {
    Unit unit(x[0]);
#ifdef SBODY_RELEASE
    auto t_start = std::chrono::steady_clock::now();
#endif
    const std::array<real_t, 9> estimate_step = {120., 60., 30., 10., 3., 1., 0.5, 0.1, 0.};
    size_t estimate_idx = 0;
    real_t this_gr_obs_time, last_gr_obs_time, gr_offset = 0.;
    real_t R = x[2] * unit.pc, inclination = x[5] * M_PI / 180., rotation = x[6] * M_PI / 180.;
    std::unique_ptr<Metric<real_t>> main_metric;
    std::unique_ptr<View<real_t>> view_ptr;
    try {
        main_metric = getMetric(metric, x[1]);
    } catch (const std::invalid_argument &e) {
        infra_context.logger_.error(e.what());
        return {py::array_t<double>(), StatusCode::INVALID_ARGUMENT};
    }
    view_ptr = std::make_unique<View<real_t>>(*main_metric, R, inclination, rotation);
    Hotspot<real_t> hotspot = initializeHotspot(unit, *main_metric, x, mode, infra_context);
    real_t t = 0., t_step = 0., t_rec = 0.1 * unit.s, t0, last_obs_time = 0., this_obs_time;
    ublas::bounded_vector<real_t, 8> last_position, this_position, cartesian_position;
    ublas::bounded_vector<real_t, 5> last_view_info, this_view_info;
    ublas::bounded_vector<real_t, 2> last_hotspot_info, this_hotspot_info;
    const real_t sin_inc = sin(inclination), cos_inc = cos(inclination);
    hotspot.position(last_position);
    if (last_position[1] < HOTSPOT_TRACE_INNER_RADIUS) {
        infra_context.logger_.error("[!] R < {}, stop!", HOTSPOT_TRACE_INNER_RADIUS);
        return {py::array_t<double>(), StatusCode::FALL_INTO_BLACK_HOLE};
    }
    real_t x0 = last_position[1] * abs(sin(last_position[2])) * cos(last_position[3]),
           z0 = last_position[1] * cosTheta(last_position[2]);
    if (gr_time_delay) {
        if (int status = view_ptr->Trace(last_position, TimeSystem::T, last_view_info, HOTSPOT_TRACE_ABS_ACC,
                                         HOTSPOT_TRACE_REL_ACC, calculate_magnification, false, infra_context);
            status != StatusCode::SUCCESS) {
            // PrintlnError("Initially trace star Error!");
            auto numpy = File(FileFormat::NUMPY,
                              fmt::format("log/hotspot/metric_{}/mode_{}/trace/{}", metric, mode,
                                          std::chrono::duration_cast<std::chrono::milliseconds>(
                                              std::chrono::system_clock::now().time_since_epoch())
                                              .count()),
                              {6});
            if (metric == 2)
                numpy.save(x[1]);
            else
                numpy.save(0.0);
            numpy.save(inclination);
            for (size_t i = 0; i < 4; ++i)
                numpy.save(last_position[i]);
            return {hotspotExit(x, infra_context), status};
        }
        t0 = last_view_info[3] / unit.s;
    }
    real_t vz0 = x0 * sin_inc + z0 * cos_inc;
    real_t h = 1.;
    int status = 0;
    size_t idx = 0, size = obs_time.size();
    auto result = py::array_t<double>(size * HOTSPOT_RETURN_WIDTH);
    double *result_ptr = result.mutable_data();
    while (true) {
#ifdef SBODY_RELEASE
        if (auto now = std::chrono::steady_clock::now(); now - t_start > std::chrono::seconds(TIME_OUT_SECONDS)) {
            // py::print("CalculateHSOrbit() timeout!");
            // auto numpy = NumPy(fmt::format("log/hotspot/metric_{}/mode_{}/timeout/{}", metric, mode,
            // chrono::duration_cast<chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()),
            // {static_cast<size_t>(x.size())}); for (int j = 0; j < x.size(); ++j) 	numpy.Save(x.at(j));
            return {hotspotExit(x, infra_context), StatusCode::COMPUTATION_TIMEOUT};
        }
#endif
        t_step += t_rec;
        if (status = hotspot.integratorApply(t, t_step, h); status != StatusCode::SUCCESS) {
            return {py::array_t<double>(), status};
        }
#ifndef SBODY_RELEASE
        hotspot.normalize();
#endif
        hotspot.position(this_position);
        if (this_position[1] < HOTSPOT_TRACE_INNER_RADIUS) {
            return {py::array_t<double>(), StatusCode::FALL_INTO_BLACK_HOLE};
        }
        this_obs_time = (t + vz0 -
                         this_position[1] * (sinTheta(this_position[2]) * cos(this_position[3]) * sin_inc +
                                             cosTheta(this_position[2]) * cos_inc)) /
                        unit.s;
        if (gr_time_delay) {
            if (this_obs_time + gr_offset + estimate_step[estimate_idx] > obs_time.at(idx)) {
                if (status = calculateHotspotOrbitWithGRDelay(
                        unit, *view_ptr, hotspot, estimate_step, estimate_idx, obs_time, idx, this_obs_time, gr_offset, t, t0,
                        t_step, t_rec, last_position, this_position, cartesian_position, last_view_info, this_view_info,
                        last_hotspot_info, this_hotspot_info, last_gr_obs_time, this_gr_obs_time, result_ptr,
                        calculate_magnification, infra_context);
                    status != StatusCode::SUCCESS)
                    return {hotspotExit(x, infra_context), status};
                if (idx >= size)
                    break;
            }
        } else if (this_obs_time > obs_time.at(idx)) {
            if (status = calculateHotspotOrbitWithRomerDelay(unit, *view_ptr, hotspot, obs_time, idx, last_obs_time,
                                                             this_obs_time, t_step, t_rec, last_position, this_position,
                                                             cartesian_position, last_view_info, result_ptr, ray_tracing,
                                                             calculate_magnification, infra_context);
                status != StatusCode::SUCCESS)
                return {hotspotExit(x, infra_context), status};
            if (idx >= size)
                break;
        }
        last_obs_time = this_obs_time;
        swap(last_position, this_position);
    }
    return {result.reshape({size, HOTSPOT_RETURN_WIDTH}), 0};
}

py::array_t<double> calculateHotspotEnergyAngularMomentumAndMinimumRadius(const std::vector<real_t> x, int metric) {
    Unit unit(x[0]);
    std::unique_ptr<Metric<real_t>> main_metric;
    InfraContext infra_context(HOTSPOT_LOG_LEVEL, HOTSPOT_STD_OUT, HOTSPOT_STD_ERR);
    auto result = py::array_t<double>(3);
    double *result_ptr = result.mutable_data();
    result_ptr[0] = result_ptr[1] = NAN;
    result_ptr[2] = 0.0;
    try {
        main_metric = getMetric(metric, x[1]);
    } catch (const std::invalid_argument &e) {
        infra_context.logger_.error(e.what());
        return py::array_t<double>();
    }
    Hotspot<real_t> hotspot = initializeHotspot(unit, *main_metric, x, 1, infra_context);
    result_ptr[0] = hotspot.energy();
    if (metric == 1) {
        result_ptr[1] = power2(hotspot.position(1)) *
                        sqrt(power2(hotspot.position(6)) + power2(std::sin(hotspot.position(2)) * hotspot.position(7)));
        real_t cubic_roots[3];
        real_t coeff = 1. / (power2(result_ptr[0]) - 1.);
        real_t L2_coeff = power2(result_ptr[1]) * coeff;
        int num_of_roots = polySolveCubic(2. * coeff, -L2_coeff, 2. * L2_coeff, cubic_roots);
        for (int i = 0; i < num_of_roots; ++i)
            if (cubic_roots[i] > 0.) {
                result_ptr[2] = cubic_roots[i];
                break;
            }
    } else
        result_ptr[1] = hotspot.angularMomentum();
    return result;
}

double hotspotFluxChi2(std::vector<real_t> &x, int metric, std::vector<real_t> &intrinsic_time,
                       std::vector<real_t> &obs_magnification, std::vector<real_t> &obs_redshift, std::vector<real_t> &obs_flux,
                       std::vector<real_t> &flux_sigma) {
    Unit unit(x.at(0));
    InfraContext infra_context(HOTSPOT_LOG_LEVEL, HOTSPOT_STD_OUT, HOTSPOT_STD_ERR);
    std::unique_ptr<Metric<real_t>> main_metric;
    try {
        main_metric = getMetric(metric, x.at(1));
    } catch (const std::invalid_argument &e) {
        infra_context.logger_.error(e.what());
        return std::numeric_limits<double>::infinity();
    }
    Hotspot<real_t> hotspot(*main_metric, TimeSystem::TAU, DynamicalSystem::LAGRANGIAN, HOTSPOT_TRACE_ABS_ACC,
                            HOTSPOT_TRACE_REL_ACC, x.at(0), x.at(1), x.at(2) * unit.s, x.at(3) * unit.s, false);
    double chi2 = 0.;
    for (size_t i = 0; i < obs_flux.size(); ++i)
        chi2 +=
            power2((fluxDensity(hotspot.spectralDensity(intrinsic_time.at(i), obs_redshift.at(i)), obs_magnification.at(i)) -
                    obs_flux.at(i)) /
                   flux_sigma.at(i));
    return chi2;
}

std::tuple<py::array_t<double>, int> hotspotDistinctChi2(const std::vector<real_t> &x, int metric, int mode, bool ray_tracing,
                                                         bool light_bending, bool gr_time_delay, bool calculate_motion_chi2,
                                                         bool calculate_flux_chi2,
                                                         const std::array<std::vector<real_t>, 7> &obs_data,
                                                         const InfraContext &infra_context) {
    auto [results, status] = calculateHotspotOrbit(x, metric, mode, ray_tracing, gr_time_delay, obs_data[OBS_TIME_IDX],
                                                   calculate_flux_chi2, infra_context);
    if (results.size() == 0)
        return {py::array_t<double>(), status};
    const size_t size = obs_data[OBS_TIME_IDX].size();
    auto result = py::array_t<double>(size);
    double *result_ptr = result.mutable_data();
    memset(result_ptr, 0, size * sizeof(double));
    if (calculate_flux_chi2)
        for (size_t i = 0; i < size; ++i)
            result_ptr[i] += power2((obs_data[OBS_FLUX_IDX].at(i) - results.at(i, 16)) / obs_data[OBS_FLUX_SIGMA_IDX].at(i));
    if (!calculate_motion_chi2)
        return {result, StatusCode::SUCCESS};
    if (light_bending)
        for (size_t i = 0; i < size; ++i) {
            result_ptr[i] +=
                power2((obs_data[OBS_RA_IDX].at(i) + results.at(i, 10) * x.at(0) * M_SUN_RG_IN_PARSEC_MAS / x.at(2) - x.at(3)) /
                       obs_data[OBS_RA_SIGMA_IDX].at(i));
            result_ptr[i] += power2(
                (obs_data[OBS_DEC_IDX].at(i) - results.at(i, 11) * x.at(0) * M_SUN_RG_IN_PARSEC_MAS / x.at(2) - x.at(4)) /
                obs_data[OBS_DEC_SIGMA_IDX].at(i));
        }
    else {
        const real_t inclination = x.at(5) * M_PI / 180., rotation = x.at(6) * M_PI / 180.;
        const real_t sin_inc = std::sin(inclination), cos_inc = std::cos(inclination);
        const real_t sin_rot = std::sin(rotation), cos_rot = std::cos(rotation);
        for (size_t i = 0; i < size; ++i) {
            const real_t inc_ra = sin_inc * results.at(i, 3) - cos_inc * results.at(i, 1);
            const real_t rot_ra = cos_rot * inc_ra - sin_rot * results.at(i, 2);
            const real_t rot_dec = sin_rot * inc_ra + cos_rot * results.at(i, 2);
            result_ptr[i] +=
                power2((obs_data[OBS_RA_IDX].at(i) + rot_ra * x.at(0) * M_SUN_RG_IN_PARSEC_MAS / x.at(2) - x.at(3)) /
                       obs_data[OBS_RA_SIGMA_IDX].at(i));
            result_ptr[i] +=
                power2((obs_data[OBS_DEC_IDX].at(i) + rot_dec * x.at(0) * M_SUN_RG_IN_PARSEC_MAS / x.at(2) - x.at(4)) /
                       obs_data[OBS_DEC_SIGMA_IDX].at(i));
        }
    }
    return {result, StatusCode::SUCCESS};
}

std::tuple<double, int> hotspotChi2(const std::vector<real_t> &x, int metric, int mode, bool ray_tracing, bool light_bending,
                                    bool gr_time_delay, double black_hole_sigma, bool calculate_motion_chi2,
                                    bool calculate_flux_chi2, const std::array<std::vector<real_t>, 7> &obs_data,
                                    const InfraContext &infra_context) {
    auto [result, status] = hotspotDistinctChi2(x, metric, mode, ray_tracing, light_bending, gr_time_delay,
                                                calculate_motion_chi2, calculate_flux_chi2, obs_data, infra_context);
    if (result.size() == 0)
        return {std::numeric_limits<double>::infinity(), status};
    double chi2 = 0.0;
    for (ssize_t i = 0; i < result.size(); ++i)
        chi2 += result.at(i);
    if (black_hole_sigma > 0.)
        chi2 += (power2(x.at(3)) + power2(x.at(4))) / power2(black_hole_sigma);
    if (isnan(chi2))
        return {std::numeric_limits<double>::infinity(), StatusCode::FAILURE};
    return {chi2, StatusCode::SUCCESS};
}

py::array_t<double> fastTrace(const real_t r_observer, const real_t theta_observer, const real_t sin_theta_observer,
                              const real_t cos_theta_observer, const real_t r_target, const real_t theta_target,
                              const real_t phi_target) {
    InfraContext infra_context(HOTSPOT_LOG_LEVEL, HOTSPOT_STD_OUT, HOTSPOT_STD_ERR);
    ublas::bounded_vector<real_t, 8> photon;
    real_t alpha = 0., beta = 0., photon_time;
    Schwarzschild<real_t>().fastTrace(r_observer, theta_observer, sin_theta_observer, cos_theta_observer, r_target,
                                      theta_target, phi_target, alpha, beta, photon, photon_time, HOTSPOT_TRACE_ABS_ACC,
                                      HOTSPOT_TRACE_REL_ACC, infra_context);
    auto res = py::array_t<double>(9);
    copy(photon.begin(), photon.end(), res.mutable_data());
    res.mutable_data()[8] = photon_time;
    return res;
}

int hotspotOrbitTest(const py::array_t<double> x, int metric, int mode) {
    Unit unit(x.at(0));
    std::vector<real_t> x_vec(x.data(), x.data() + x.size());
    InfraContext infra_context(HOTSPOT_LOG_LEVEL, HOTSPOT_STD_OUT, HOTSPOT_STD_ERR);
    std::unique_ptr<Metric<real_t>> main_metric;
    std::unique_ptr<View<real_t>> view_ptr;
    try {
        main_metric = getMetric(metric, x.at(1));
    } catch (const std::invalid_argument &e) {
        infra_context.logger_.error(e.what());
        return StatusCode::FAILURE;
    }
    int status = StatusCode::SUCCESS;
    Hotspot<real_t> hotspot = initializeHotspot(unit, *main_metric, x_vec, mode, infra_context);
    real_t t = 0., t_step = 0., t_rec = 10. * unit.s;
    real_t h = 1.0;
    for (int i = 0; i < 100; ++i) {
        t_step += t_rec;
        if (status = hotspot.integratorApply(t, t_step, h); status != StatusCode::SUCCESS) {
            infra_context.logger_.error("[!] IntegratorApply status =", status);
            return StatusCode::FAILURE;
        }
        if (hotspot.position(1) < HOTSPOT_TRACE_INNER_RADIUS)
            break;
#ifndef SBODY_RELEASE
        hotspot.normalize();
#endif
        infra_context.logger_.info("r:", hotspot.position(1), "Energy:", hotspot.energy(),
                                   "Angular Momentum:", hotspot.angularMomentum());
    }
    return StatusCode::SUCCESS;
}

class HotspotEvent {
  protected:
    const int metric_;
    const int motion_mode_;
    const bool ray_tracing_;
    const bool light_bending_;
    const bool gr_time_delay_;
    const bool black_hole_constraint_;
    const size_t value_dim_;
    const std::vector<real_t> fix_value_;
    const std::vector<uint8_t> fix_mask_;
    const std::vector<real_t> lower_bound_;
    const std::vector<real_t> upper_bound_;
    std::vector<real_t> embed_x_;
    const size_t data_num_;
    const std::array<std::vector<real_t>, 7> obs_data_;
    // const std::vector<real_t> obs_time_;
    // const std::vector<real_t> flux_;
    // const std::vector<real_t> flux_sigma_;
    // const std::vector<real_t> ra_;
    // const std::vector<real_t> ra_sigma_;
    // const std::vector<real_t> dec_;
    // const std::vector<real_t> dec_sigma_;
    const double black_hole_sigma_;
    const InfraContext infra_context_;

  public:
    HotspotEvent(int metric, int mode, int gr_switch, std::vector<real_t> fix_value, std::vector<uint8_t> fix_mask,
                 std::vector<real_t> lower_bound, std::vector<real_t> upper_bound, std::array<std::vector<real_t>, 7> obs_data,
                 double black_hole_sigma)
        : metric_(metric), motion_mode_(mode), ray_tracing_(gr_switch > 0), light_bending_((gr_switch & 2) > 0),
          gr_time_delay_((gr_switch & 4) > 0), black_hole_constraint_((gr_switch & 16) > 0), value_dim_(fix_value.size()),
          fix_value_(std::move(fix_value)), fix_mask_(std::move(fix_mask)), lower_bound_(std::move(lower_bound)),
          upper_bound_(std::move(upper_bound)), embed_x_(value_dim_), data_num_(obs_data[0].size()),
          obs_data_(std::move(obs_data)), black_hole_sigma_(black_hole_sigma),
          infra_context_(HOTSPOT_LOG_LEVEL, HOTSPOT_STD_OUT, HOTSPOT_STD_ERR) {
        assert(value_dim_ == fix_mask_.size());
        assert(value_dim_ == lower_bound_.size());
        assert(value_dim_ == upper_bound_.size());
        for (size_t i = 0; i < value_dim_; ++i) {
            if (!fix_mask_[i])
                continue;
            assert(fix_value_[i] >= lower_bound_[i]);
            assert(fix_value_[i] < upper_bound_[i]);
            embed_x_[i] = fix_value_[i];
        }
        for (size_t i = 1; i < 7; ++i)
            assert(data_num_ == obs_data_[i].size());
    };

    bool embedCheck(const py::array_t<double> x) {
        auto x_uncheck = x.unchecked<1>();
        ssize_t x_idx = 0;
        for (size_t i = 0; i < value_dim_; ++i)
            if (!fix_mask_[i]) {
                embed_x_[i] = x_uncheck(x_idx++);
                if (embed_x_[i] < lower_bound_[i] || embed_x_[i] >= upper_bound_[i])
                    return true;
            }
        return false;
    }

    py::array_t<double> calculateFullOrbit(const py::array_t<double> x, double t1) {
        if (embedCheck(x))
            return py::array_t<double>();
        return calculateHotspotFullOrbit(embed_x_, metric_, motion_mode_, gr_time_delay_, t1);
    }
    std::tuple<py::array_t<double>, int> calculateOrbit(const py::array_t<double> x, bool calculate_magnification) {
        if (embedCheck(x))
            return {py::array_t<double>(), StatusCode::INVALID_ARGUMENT};
        return calculateHotspotOrbit(embed_x_, metric_, motion_mode_, ray_tracing_, gr_time_delay_, obs_data_[0],
                                     calculate_magnification, infra_context_);
    }
    double motionChi2FromOrbitResult(const py::array_t<double> orbit_result) {
        if (orbit_result.size() == 0)
            return std::numeric_limits<double>::infinity();
        double chi2 = 0.0;
        if (light_bending_)
            for (size_t i = 0; i < data_num_; ++i) {
                chi2 += power2((obs_data_[OBS_RA_IDX].at(i) +
                                orbit_result.at(i, 10) * embed_x_.at(0) * M_SUN_RG_IN_PARSEC_MAS / embed_x_.at(2) -
                                embed_x_.at(3)) /
                               obs_data_[OBS_RA_SIGMA_IDX].at(i));
                chi2 += power2((obs_data_[OBS_DEC_IDX].at(i) -
                                orbit_result.at(i, 11) * embed_x_.at(0) * M_SUN_RG_IN_PARSEC_MAS / embed_x_.at(2) -
                                embed_x_.at(4)) /
                               obs_data_[OBS_DEC_SIGMA_IDX].at(i));
            }
        else {
            infra_context_.logger_.error("Inclination and rotation is ignored!");
            real_t inclination = embed_x_[5] * M_PI / 180., rotation = embed_x_[6] * M_PI / 180.;
            real_t sin_inc = std::sin(inclination), cos_inc = std::cos(inclination);
            real_t sin_rot = std::sin(rotation), cos_rot = std::cos(rotation);
            for (size_t i = 0; i < data_num_; ++i) {
                const real_t inc_ra = sin_inc * orbit_result.at(i, 3) - cos_inc * orbit_result.at(i, 1);
                const real_t rot_ra = cos_rot * inc_ra - sin_rot * orbit_result.at(i, 2);
                const real_t rot_dec = sin_rot * inc_ra + cos_rot * orbit_result.at(i, 2);
                chi2 += power2((obs_data_[OBS_RA_IDX].at(i) +
                                rot_ra * embed_x_.at(0) * M_SUN_RG_IN_PARSEC_MAS / embed_x_.at(2) - embed_x_.at(3)) /
                               obs_data_[OBS_RA_SIGMA_IDX].at(i));
                chi2 += power2((obs_data_[OBS_DEC_IDX].at(i) +
                                rot_dec * embed_x_.at(0) * M_SUN_RG_IN_PARSEC_MAS / embed_x_.at(2) - embed_x_.at(4)) /
                               obs_data_[OBS_DEC_SIGMA_IDX].at(i));
            }
        }
        return chi2;
    }
    std::tuple<double, int> motionLogProb(const py::array_t<double> x) {
        if (embedCheck(x))
            return {-std::numeric_limits<double>::infinity(), StatusCode::INVALID_ARGUMENT};
        auto [chi2, status] =
            hotspotChi2(embed_x_, metric_, motion_mode_, ray_tracing_, light_bending_, gr_time_delay_,
                        black_hole_constraint_ ? black_hole_sigma_ : 0.0, true, false, obs_data_, infra_context_);
        if (status != StatusCode::SUCCESS)
            return {-std::numeric_limits<double>::infinity(), status};
        if (isnan(chi2))
            return {-std::numeric_limits<double>::infinity(), StatusCode::FAILURE};
        return {-0.5 * chi2, StatusCode::SUCCESS};
    }
    std::tuple<double, int> fluxLogProb(const py::array_t<double> x) {
        if (embedCheck(x))
            return {-std::numeric_limits<double>::infinity(), StatusCode::INVALID_ARGUMENT};
        auto [chi2, status] =
            hotspotChi2(embed_x_, metric_, motion_mode_, ray_tracing_, light_bending_, gr_time_delay_,
                        black_hole_constraint_ ? black_hole_sigma_ : 0.0, false, true, obs_data_, infra_context_);
        if (status != StatusCode::SUCCESS)
            return {-std::numeric_limits<double>::infinity(), status};
        if (isnan(chi2))
            return {-std::numeric_limits<double>::infinity(), StatusCode::FAILURE};
        return {-0.5 * chi2, StatusCode::SUCCESS};
    }
    std::tuple<double, int> totalLogProb(const py::array_t<double> x) {
        if (embedCheck(x))
            return {-std::numeric_limits<double>::infinity(), StatusCode::INVALID_ARGUMENT};
        auto [chi2, status] =
            hotspotChi2(embed_x_, metric_, motion_mode_, ray_tracing_, light_bending_, gr_time_delay_,
                        black_hole_constraint_ ? black_hole_sigma_ : 0.0, true, true, obs_data_, infra_context_);
        if (status != StatusCode::SUCCESS)
            return {-std::numeric_limits<double>::infinity(), status};
        if (isnan(chi2))
            return {-std::numeric_limits<double>::infinity(), StatusCode::FAILURE};
        return {-0.5 * chi2, StatusCode::SUCCESS};
    }
};

std::vector<size_t> getEventDataNumList(const std::vector<std::array<std::vector<real_t>, 7>> &obs_data_list) {
    std::vector<size_t> data_num_list;
    for (const auto &obs_data : obs_data_list)
        data_num_list.push_back(obs_data[0].size());
    return data_num_list;
}

class HotspotCombinationEvent {
  protected:
    const int metric_;
    const int motion_mode_;
    const bool ray_tracing_;
    const bool light_bending_;
    const bool gr_time_delay_;
    const bool black_hole_constraint_;
    // const int gr_switch_;
    const size_t value_dim_;
    const std::vector<real_t> fix_value_;
    const std::vector<uint8_t> fix_mask_;
    const std::vector<real_t> lower_bound_;
    const std::vector<real_t> upper_bound_;
    const size_t event_num_;
    std::vector<std::vector<real_t>> embed_x_list_;
    const size_t first_event_dim_;
    const size_t additional_dim_;
    const std::vector<uint8_t> combination_mask_;
    const std::vector<size_t> data_num_list_;
    const std::vector<std::array<std::vector<real_t>, 7>> obs_data_list_;
    const double black_hole_sigma_;
    const InfraContext infra_context_;

  public:
    HotspotCombinationEvent(int metric, int mode, int gr_switch, std::vector<real_t> fix_value, std::vector<uint8_t> fix_mask,
                            std::vector<real_t> lower_bound, std::vector<real_t> upper_bound,
                            std::vector<uint8_t> combination_mask,
                            std::vector<std::array<std::vector<real_t>, 7>> obs_data_list, double black_hole_sigma)
        : metric_(metric), motion_mode_(mode), ray_tracing_(gr_switch > 0), light_bending_((gr_switch & 2) > 0),
          gr_time_delay_((gr_switch & 4) > 0), black_hole_constraint_((gr_switch & 16) > 0), value_dim_(fix_value.size()),
          fix_value_(std::move(fix_value)), fix_mask_(std::move(fix_mask)), lower_bound_(std::move(lower_bound)),
          upper_bound_(std::move(upper_bound)), event_num_(obs_data_list.size()),
          embed_x_list_(event_num_, std::vector<real_t>(value_dim_)), first_event_dim_(combination_mask.size()),
          additional_dim_(std::accumulate(combination_mask.begin(), combination_mask.end(), 0)),
          combination_mask_(std::move(combination_mask)), data_num_list_(getEventDataNumList(obs_data_list)),
          obs_data_list_(std::move(obs_data_list)), black_hole_sigma_(black_hole_sigma),
          infra_context_(HOTSPOT_LOG_LEVEL, HOTSPOT_STD_OUT, HOTSPOT_STD_ERR) {
        assert(value_dim_ == fix_mask_.size());
        assert(value_dim_ == lower_bound_.size());
        assert(value_dim_ == upper_bound_.size());
        for (size_t i = 0; i < value_dim_; ++i) {
            if (!fix_mask_[i])
                continue;
            assert(fix_value_[i] >= lower_bound_[i]);
            assert(fix_value_[i] < upper_bound_[i]);
            for (size_t j = 0; j < event_num_; ++j)
                embed_x_list_[j][i] = fix_value_[i];
        }
        for (size_t i = 0; i < event_num_; ++i)
            for (size_t j = 1; j < 7; ++j)
                assert(data_num_list_[i] == obs_data_list_[i][j].size());
    };

    bool embedCheck(const py::array_t<double> x) {
        ssize_t x_idx = 0, embed_idx = first_event_dim_ - additional_dim_;
        auto x_uncheck = x.unchecked<1>();
        for (size_t i = 0; i < value_dim_; ++i)
            if (!fix_mask_[i]) {
                embed_x_list_[0][i] = x_uncheck(x_idx);
                if (embed_x_list_[0][i] < lower_bound_[i] || embed_x_list_[0][i] >= upper_bound_[i])
                    return true;
                if (combination_mask_[x_idx]) {
                    for (size_t j = 1; j < event_num_; ++j)
                        embed_x_list_[j][i] = embed_x_list_[0][i];
                } else {
                    for (size_t j = 1; j < event_num_; ++j) {
                        embed_x_list_[j][i] = x_uncheck(additional_dim_ * j + embed_idx);
                        if (embed_x_list_[j][i] < lower_bound_[i] || embed_x_list_[j][i] >= upper_bound_[i])
                            return true;
                    }
                    ++embed_idx;
                }
                ++x_idx;
            }
        return false;
    }

    bool conicalFlip(int conical_switch) {
        if (conical_switch == 0)
            return false;
        for (size_t i = 1; i < event_num_; ++i) {
            if (conical_switch & 1)
                embed_x_list_[i][8] = 180. - embed_x_list_[i][8];
            if (conical_switch >>= 1; conical_switch == 0)
                return false;
        }
        return true;
    }

    void printEmbed() {
        for (size_t i = 0; i < event_num_; ++i) {
            py::print("Event", i);
            for (size_t j = 0; j < value_dim_; ++j)
                py::print("  ", embed_x_list_[i][j]);
        }
    }

    py::array_t<double> calculateFullOrbit(const py::array_t<double> x, size_t event_idx, double t1) {
        if (embedCheck(x))
            return py::array_t<double>();
        return calculateHotspotFullOrbit(embed_x_list_[event_idx], metric_, motion_mode_, gr_time_delay_, t1);
    }

    std::tuple<py::array_t<double>, int> calculateOrbit(const py::array_t<double> x, size_t event_idx,
                                                        bool calculate_magnification) {
        if (embedCheck(x))
            return {py::array_t<double>(), StatusCode::INVALID_ARGUMENT};
        return calculateHotspotOrbit(embed_x_list_[event_idx], metric_, motion_mode_, ray_tracing_, gr_time_delay_,
                                     obs_data_list_[event_idx][OBS_TIME_IDX], calculate_magnification, infra_context_);
    }

    std::tuple<double, int> motionLogProb(const py::array_t<double> x) {
        if (embedCheck(x))
            return {-std::numeric_limits<double>::infinity(), StatusCode::INVALID_ARGUMENT};
        double chi2 = 0.0;
        for (size_t i = 0; i < event_num_; ++i) {
            auto [chi2_i, status] = hotspotChi2(embed_x_list_[i], metric_, motion_mode_, ray_tracing_, light_bending_,
                                                gr_time_delay_, 0.0, true, false, obs_data_list_[i], infra_context_);
            if (status != StatusCode::SUCCESS)
                return {-std::numeric_limits<double>::infinity(), status};
            chi2 += chi2_i;
        }
        if (black_hole_constraint_)
            chi2 += (power2(embed_x_list_[0].at(3)) + power2(embed_x_list_[0].at(4))) / power2(black_hole_sigma_);
        return {-0.5 * chi2, StatusCode::SUCCESS};
    }
    std::tuple<double, int> fluxLogProb(const py::array_t<double> x) {
        if (embedCheck(x))
            return {-std::numeric_limits<double>::infinity(), StatusCode::INVALID_ARGUMENT};
        double chi2 = 0.0;
        for (size_t i = 0; i < event_num_; ++i) {
            auto [chi2_i, status] = hotspotChi2(embed_x_list_[i], metric_, motion_mode_, ray_tracing_, light_bending_,
                                                gr_time_delay_, 0.0, false, true, obs_data_list_[i], infra_context_);
            if (status != StatusCode::SUCCESS)
                return {-std::numeric_limits<double>::infinity(), status};
            chi2 += chi2_i;
        }
        return {-0.5 * chi2, StatusCode::SUCCESS};
    }
    std::tuple<double, int> totalLogProb(const py::array_t<double> x) {
        if (embedCheck(x))
            return {-std::numeric_limits<double>::infinity(), StatusCode::INVALID_ARGUMENT};
        double chi2 = 0.0;
        for (size_t i = 0; i < event_num_; ++i) {
            auto [chi2_i, status] = hotspotChi2(embed_x_list_[i], metric_, motion_mode_, ray_tracing_, light_bending_,
                                                gr_time_delay_, 0.0, true, true, obs_data_list_[i], infra_context_);
            if (status != StatusCode::SUCCESS)
                return {-std::numeric_limits<double>::infinity(), status};
            chi2 += chi2_i;
        }
        if (black_hole_constraint_)
            chi2 += (power2(embed_x_list_[0].at(3)) + power2(embed_x_list_[0].at(4))) / power2(black_hole_sigma_);
        return {-0.5 * chi2, StatusCode::SUCCESS};
    }
};

void bindHotspotModule(py::module_ &m) {
    py::class_<HotspotEvent>(m, "HotspotEvent")
        .def(py::init<int, int, int, std::vector<real_t>, std::vector<uint8_t>, std::vector<real_t>, std::vector<real_t>,
                      std::array<std::vector<real_t>, 7>, double>())
        .def("calculate_full_orbit", &HotspotEvent::calculateFullOrbit, py::arg("x"), py::arg("t1"))
        .def("calculate_orbit", &HotspotEvent::calculateOrbit, py::arg("x"), py::arg("calculate_magnification"))
        .def("motion_log_prob", &HotspotEvent::motionLogProb, py::arg("x"))
        // .def("motionChi2FromOrbitResult", &HotspotEvent::motionChi2FromOrbitResult, py::arg("orbit_result"))
        .def("flux_log_prob", &HotspotEvent::fluxLogProb, py::arg("x"))
        .def("total_log_prob", &HotspotEvent::totalLogProb, py::arg("x"));

    py::class_<HotspotCombinationEvent>(m, "HotspotCombinationEvent")
        .def(py::init<int, int, int, std::vector<real_t>, std::vector<uint8_t>, std::vector<real_t>, std::vector<real_t>,
                      std::vector<uint8_t>, std::vector<std::array<std::vector<real_t>, 7>>, double>())
        .def("embed_check", &HotspotCombinationEvent::embedCheck, py::arg("x"))
        .def("print_embed", &HotspotCombinationEvent::printEmbed)
        .def("calculate_full_orbit", &HotspotCombinationEvent::calculateFullOrbit, py::arg("x"), py::arg("event_idx"),
             py::arg("t1"))
        .def("calculate_orbit", &HotspotCombinationEvent::calculateOrbit, py::arg("x"), py::arg("event_idx"),
             py::arg("calculate_magnification"))
        .def("motion_log_prob", &HotspotCombinationEvent::motionLogProb, py::arg("x"))
        .def("flux_log_prob", &HotspotCombinationEvent::fluxLogProb, py::arg("x"))
        .def("total_log_prob", &HotspotCombinationEvent::totalLogProb, py::arg("x"));

    m.def("calculate_hotspot_full_orbit", &calculateHotspotFullOrbit, R"pbdoc(
    	Get trajectory, redshift, apparent position with orbital parameters
    )pbdoc");

    m.def("calculate_hotspot_orbit", &calculateHotspotOrbit, R"pbdoc(
		Get trajectory, redshift, apparent position with orbital parameters
	)pbdoc");

    m.def("hotspot_chi2", &hotspotChi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

    m.def("hotspot_distinct_chi2", &hotspotDistinctChi2, R"pbdoc(
		Get MCMC log probability
	)pbdoc");

    m.def("hotspot_flux_chi2", &hotspotFluxChi2, R"pbdoc(
        Get MCMC log probability with only flux data
    )pbdoc");

    m.def("calculate_hotspot_energy_angular_momentum_and_minimum_radius",
          &calculateHotspotEnergyAngularMomentumAndMinimumRadius,
          R"pbdoc(
		Calculate Hotspot Energy Angular Momentum And Minimum Radius for geodesic orbits.
	)pbdoc");

    m.def("fast_trace", &fastTrace, R"pbdoc(
		Schwarzschild fast trace
	)pbdoc");
    m.def("hotspot_orbit_test", &hotspotOrbitTest, R"pbdoc(
		HotspotOrbitTest
	)pbdoc");
}
