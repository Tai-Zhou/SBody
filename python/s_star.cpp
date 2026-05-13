/**
 * @file s_star.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <cassert>
#include <csignal>
#include <memory>

#include <boost/algorithm/algorithm.hpp>
#include <fmt/core.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include "infra/infra_context.hpp"
#include "infra/log.hpp"
#include "infra/status.hpp"
#include "math/geometry.hpp"
#include "math/interpolation.hpp"
#include "metric/post_newtonian.hpp"
#include "metric/schwarzschild.hpp"
#include "object/particle.hpp"
#include "units.hpp"
#include "view/view.hpp"

namespace py = pybind11;
namespace ublas = boost::numeric::ublas;
using namespace std;
using namespace sbody;

constexpr double M_SUN_RG_IN_PARSEC_MAS = 1.0 / (20896825544594.496 * mas_in_rad);

constexpr double S_STAR_TRACE_ABS_ACC = 1e-12;
constexpr double S_STAR_TRACE_REL_ACC = 1e-12;

double CalculatePericenterTime(double mass, int metric, double fSP, double R, double t_apo, double a, double e,
                               double inclination, double ascending_node, double periapsis) {
    Unit unit(mass);
    R *= unit.pc;
    a *= R * mas_in_rad;
    inclination *= M_PI / 180.;
    ascending_node *= M_PI / 180.;
    periapsis *= M_PI / 180.;
    double t = (t_apo - 2002.) * unit.yr, tStep = 16.3 * unit.yr, tRec = 0.000001 * unit.yr;
    unique_ptr<Metric<double>> main_metric;
    unique_ptr<View<double>> view_ptr;
    InfraContext infra_context;
    if (metric == 0)
        main_metric = make_unique<PN1<double>>(fSP);
    else {
        main_metric = make_unique<Schwarzschild<double>>();
        view_ptr = make_unique<View<double>>(*main_metric, R, 0., 0.);
    }
    auto star_0 =
        Particle(*main_metric, TimeSystem::T, DynamicalSystem::LAGRANGIAN, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC, false);
    if (metric == 2)
        star_0.initializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
    else
        star_0.initializeKeplerian(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
    double h = 1., last_radius = std::numeric_limits<double>::infinity(), t0;
    ublas::bounded_vector<double, 8> position;
    ublas::bounded_vector<double, 5> record;
    star_0.position(position);
    if (metric == 0)
        t0 = position[1] * cosTheta(position[2]);
    else {
        view_ptr->Trace(position, TimeSystem::T, record, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC, false, false,
                        infra_context);
        t0 = record[3] / unit.s;
    }
    int status = 0;
    while (true) {
        status = star_0.integratorApply(t, tStep, h);
        if (status > 0)
            infra_context.logger_.error("main status = {}", status);
        star_0.position(position);
        if (position[1] > last_radius) {
            if (metric == 0)
                return (tStep + t0 - position[1] * cosTheta(position[2])) / unit.yr + 2002.;
            else {
                view_ptr->Trace(position, TimeSystem::T, record, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC, false, true,
                                infra_context);
                return (tStep + (t0 * unit.s - record[3])) / unit.yr + 2002.;
            }
        }
        last_radius = position[1];
        tStep += tRec;
    }
}

py::array_t<double> CalculateFullStarOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e,
                                           double inclination, double ascending_node, double periapsis, double t1) {
    Unit unit(mass);
    R *= unit.pc;
    a *= R * mas_in_rad;
    inclination *= M_PI / 180.;
    ascending_node *= M_PI / 180.;
    periapsis *= M_PI / 180.;
    size_t tStepNumber = 1000;
    double t = (tp - 2002.) * unit.yr - M_PI * sqrt(power3(a)), tStep = 0., tRec = (t1 - 2002.) * unit.yr / tStepNumber;
    unique_ptr<Metric<double>> main_metric;
    unique_ptr<View<double>> view_ptr;
    InfraContext infra_context;
    if (metric == 0) {
        main_metric = make_unique<PN1<double>>(fSP);
    } else if (metric == 3) {
        main_metric = make_unique<PostNewtonian<double>>(3);
    } else {
        main_metric = make_unique<Schwarzschild<double>>();
        view_ptr = make_unique<View<double>>(*main_metric, R, 0., 0.);
    }
    TimeSystem star_time = TimeSystem::T;
    Particle<double> star_0(*main_metric, star_time, DynamicalSystem::LAGRANGIAN, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC,
                            false);
    if (metric == 2)
        star_0.initializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
    else
        star_0.initializeKeplerian(a, e, inclination, periapsis, ascending_node, M_PI, 0., 0.);
    double h = -1.;
    int status = 0;
    if (status = star_0.integratorApply(t, 0., h); status != 0)
        py::print("[!] IntegratorApply status =", status);
    h = 1.;
    auto result = py::array_t<double>(tStepNumber * 14);
    double *result_ptr = result.mutable_data();
    ublas::bounded_vector<double, 8> position;
    ublas::bounded_vector<double, 5> record;
    for (size_t i = 0; i < tStepNumber; ++i) {
        tStep += tRec;
        status = star_0.integratorApply(t, tStep, h);
        if (status > 0)
            infra_context.logger_.error("main status = {}", status);
        star_0.position(position);
        if (metric)
            view_ptr->Trace(position, star_time, record, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC, false, true,
                            infra_context);
        sphericalToCartesian(position);
        copy(position.begin(), position.end(), result_ptr);
        record[3] /= unit.s;
        copy(record.begin(), record.end(), result_ptr + 10);
        result_ptr[8] = t / unit.s;
        const double delta_epsilon = 1. - 2. / hypot(result_ptr[1], result_ptr[2], result_ptr[3]);
        result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) /
                             sqrt(delta_epsilon - hypotSquare(result_ptr[5], result_ptr[6], result_ptr[7])) -
                         1.) *
                        299792.458;
        result_ptr += 14;
    }
    return result.reshape({tStepNumber, 14UL});
}

py::array_t<double> CalculateStarOrbit(double mass, int metric, double fSP, double R, double tp, double a, double e,
                                       double inclination, double ascending_node, double periapsis, bool ray_tracing,
                                       bool gr_time_delay, const py::array_t<double> obs_time) {
    Unit unit(mass);
    R *= unit.pc;
    a *= R * mas_in_rad;
    inclination *= M_PI / 180.;
    ascending_node *= M_PI / 180.;
    periapsis *= M_PI / 180.;
    double t = (tp - 2002.) * unit.yr - M_PI * sqrt(power3(a)), tStep = 0., tRec = 0.001 * unit.yr;
    unique_ptr<Metric<double>> main_metric;
    unique_ptr<View<double>> view_ptr;
    InfraContext infra_context;
    if (metric == 0) {
        main_metric = make_unique<PN1<double>>(fSP);
        ray_tracing = false;
        gr_time_delay = false;
    } else if (metric == 3) {
        main_metric = make_unique<PostNewtonian<double>>(3);
    } else {
        main_metric = make_unique<Schwarzschild<double>>();
        view_ptr = make_unique<View<double>>(*main_metric, R, M_PI_4, 0.);
    }
    TimeSystem star_time = TimeSystem::T;
    Particle star_0(*main_metric, star_time, DynamicalSystem::LAGRANGIAN, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC, false);
    double z0, t0, last_obs_time = 2002., this_obs_time;
    ublas::bounded_vector<double, 8> last_position, this_position, cartesian_position;
    ublas::bounded_vector<double, 5> last_view_info, this_view_info;
    if (metric == 2)
        star_0.initializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, 0., M_PI_4, 0.);
    else
        star_0.initializeKeplerian(a, e, inclination, periapsis, ascending_node, 0., M_PI_4, 0.);
    star_0.position(last_position);
    z0 = last_position[1] * cosTheta(last_position[2]);
    if (gr_time_delay) {
        if (view_ptr->Trace(last_position, star_time, last_view_info, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC, false, true,
                            infra_context) != StatusCode::SUCCESS) {
            infra_context.logger_.error("Initially trace star Error!");
            return py::array_t<double>();
        }
        t0 = last_view_info[3] / unit.s;
    }
    if (metric == 2)
        star_0.initializeKeplerianHarmonic(a, e, inclination, periapsis, ascending_node, M_PI, M_PI_4, 0.);
    else
        star_0.initializeKeplerian(a, e, inclination, periapsis, ascending_node, M_PI, M_PI_4, 0.);
    double h = -1.;
    int status = 0;
    if (status = star_0.integratorApply(t, 0., h); status != 0)
        py::print("[!] IntegratorApply status =", status);
    h = 1.;
    size_t idx = 0, size = obs_time.size();
    auto result = py::array_t<double>(size * 14);
    double *result_ptr = result.mutable_data();
    while (true) {
        tStep += tRec;
        if (status = star_0.integratorApply(t, tStep, h); status != 0)
            py::print("[!] IntegratorApply status =", status);
        star_0.position(this_position);
        if (metric == 2)
            this_obs_time = (t + z0 - (this_position[1] - 1.) * cosTheta(this_position[2])) / unit.yr + 2002.;
        else
            this_obs_time = (t + z0 - this_position[1] * cosTheta(this_position[2])) / unit.yr + 2002.;
        if (this_obs_time > obs_time.at(idx)) {
            if (gr_time_delay) {
                if (view_ptr->Trace(last_position, star_time, last_view_info, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC, false,
                                    true, infra_context) != StatusCode::SUCCESS) {
                    infra_context.logger_.error("Trace last_position Error!");
                    return py::array_t<double>();
                }
                if (view_ptr->Trace(this_position, star_time, this_view_info, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC, false,
                                    true, infra_context) != StatusCode::SUCCESS) {
                    infra_context.logger_.error("Trace this_position Error!");
                    return py::array_t<double>();
                }
                double this_gr_obs_time = (t + (t0 * unit.s - this_view_info[3])) / unit.yr + 2002.;
                double last_gr_obs_time = (t - tRec + (t0 * unit.s - last_view_info[3])) / unit.yr + 2002.;
                interpolateSphericalPositionToCartesian(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, last_position,
                                                        this_position, cartesian_position);
                copy(cartesian_position.begin(), cartesian_position.end(), result_ptr);
                result_ptr[8] =
                    linearInterpolation(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, tStep - tRec, tStep) / unit.s;
                linearInterpolation(obs_time.at(idx), last_gr_obs_time, this_gr_obs_time, last_view_info, this_view_info,
                                    result_ptr + 10);
            } else {
                interpolateSphericalPositionToCartesian(obs_time.at(idx), last_obs_time, this_obs_time, last_position,
                                                        this_position, cartesian_position);
                if (metric != 2)
                    copy(cartesian_position.begin(), cartesian_position.end(), result_ptr);
                if (ray_tracing) {
                    cartesianToSpherical(cartesian_position);
                    last_position = cartesian_position;
                    if (view_ptr->Trace(last_position, star_time, last_view_info, S_STAR_TRACE_ABS_ACC, S_STAR_TRACE_REL_ACC,
                                        false, true, infra_context) != StatusCode::SUCCESS) {
                        infra_context.logger_.error("Trace star position Error!");
                        return py::array_t<double>();
                    }
                    last_view_info[3] /= unit.s;
                    copy(last_view_info.begin(), last_view_info.end(), result_ptr + 10);
                }
                if (metric == 2) {
                    cartesianToSpherical(cartesian_position);
                    cartesian_position[1] -= 1.;
                    sphericalToCartesian(cartesian_position);
                    copy(cartesian_position.begin(), cartesian_position.end(), result_ptr);
                }
                result_ptr[8] =
                    linearInterpolation(obs_time.at(idx), last_obs_time, this_obs_time, tStep - tRec, tStep) / unit.s;
            }
            const double delta_epsilon = 1. - 2. / hypot(result_ptr[1], result_ptr[2], result_ptr[3]);
            result_ptr[9] = ((1. - result_ptr[7] / sqrt(delta_epsilon)) /
                                 sqrt(delta_epsilon - hypotSquare(result_ptr[5], result_ptr[6], result_ptr[7])) -
                             1.) *
                            299792.458;
            result_ptr += 14;
            if (++idx >= size)
                break;
        }
        last_obs_time = this_obs_time;
        swap(last_position, this_position);
        swap(last_view_info, this_view_info);
    }
    return result.reshape({size, 14UL});
}

double StarChi2(py::array_t<double> x, int metric, int gr_switch, py::array_t<double> obs_time,
                py::array_t<double> obs_redshift, py::array_t<double> redshift_sigma, py::array_t<double> obs_ra,
                py::array_t<double> ra_sigma, py::array_t<double> obs_dec, py::array_t<double> dec_sigma) {
    if (x.at(9) > 0.99)
        return std::numeric_limits<double>::infinity();
    if (metric == 0 && gr_switch > 0) {
        x.mutable_at(7) = 1.;
        gr_switch = 0;
    }
    auto obs_data = CalculateStarOrbit(x.at(0), metric, x.at(7), x.at(3), x.at(13), x.at(8), x.at(9), x.at(10), x.at(11),
                                       x.at(12), gr_switch > 0, (gr_switch & 4) > 0, obs_time);
    double redshift_prob = 0., ra_prob = 0., dec_prob = 0.;
    const int size = obs_redshift.size();
    if (gr_switch & 1) {
        for (int i = 0; i < size; ++i)
            if (obs_redshift.at(i) != std::numeric_limits<double>::infinity())
                redshift_prob +=
                    power2((obs_redshift.at(i) - (obs_data.at(i, 12) - 1.) * 299792.458 - x.at(6)) / redshift_sigma.at(i));
    } else
        for (int i = 0; i < size; ++i)
            if (obs_redshift.at(i) != std::numeric_limits<double>::infinity())
                redshift_prob += power2((obs_redshift.at(i) - obs_data.at(i, 9) - x.at(6)) / redshift_sigma.at(i));
    if (gr_switch & 2)
        for (int i = 0; i < size; ++i) {
            if (obs_ra.at(i) != std::numeric_limits<double>::infinity())
                ra_prob += power2((obs_ra.at(i) + obs_data.at(i, 10) * x.at(0) * M_SUN_RG_IN_PARSEC_MAS / x.at(3) - x.at(2) -
                                   (obs_time.at(i) - 2009.02) * x.at(5)) /
                                  ra_sigma.at(i));
            if (obs_dec.at(i) != std::numeric_limits<double>::infinity())
                dec_prob += power2((obs_dec.at(i) - obs_data.at(i, 11) * x.at(0) * M_SUN_RG_IN_PARSEC_MAS / x.at(3) - x.at(1) -
                                    (obs_time.at(i) - 2009.02) * x.at(4)) /
                                   dec_sigma.at(i));
        }
    else
        for (int i = 0; i < size; ++i) {
            if (obs_ra.at(i) != std::numeric_limits<double>::infinity())
                ra_prob += power2((obs_ra.at(i) + obs_data.at(i, 2) * x.at(0) * M_SUN_RG_IN_PARSEC_MAS / x.at(3) - x.at(2) -
                                   (obs_time.at(i) - 2009.02) * x.at(5)) /
                                  ra_sigma.at(i));
            if (obs_dec.at(i) != std::numeric_limits<double>::infinity())
                dec_prob += power2((obs_dec.at(i) + obs_data.at(i, 1) * x.at(0) * M_SUN_RG_IN_PARSEC_MAS / x.at(3) - x.at(1) -
                                    (obs_time.at(i) - 2009.02) * x.at(4)) /
                                   dec_sigma.at(i));
        }
    double chi2 = redshift_prob + ra_prob + dec_prob;
    if (isnan(chi2)) {
        if (isnan(redshift_prob))
            py::print("[!] redshift_prob is NaN");
        if (isnan(ra_prob))
            py::print("[!] ra_prob is NaN");
        if (isnan(dec_prob))
            py::print("[!] dec_prob is NaN");
        auto obs_data_test = CalculateStarOrbit(x.at(0), metric, x.at(7), x.at(3), x.at(13), x.at(8), x.at(9), x.at(10),
                                                x.at(11), x.at(12), gr_switch > 0, (gr_switch & 4) > 0, obs_time);
        for (int i = 0; i < 14; ++i)
            py::print(x.at(i));
    }
    return chi2;
}

void bindSStarModule(py::module_ &m) {
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
}
