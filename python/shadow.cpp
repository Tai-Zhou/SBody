/**
 * @file shadow.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-04-02
 *
 * @copyright Copyright (c) 2026
 */

#include <cassert>
#include <csignal>
#include <vector>

#include <boost/algorithm/algorithm.hpp>
#include <fmt/core.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include "infra/progress.hpp"
#include "metric/kerr.hpp"
#include "units.hpp"

namespace py = pybind11;
using namespace std;
using namespace sbody;

constexpr double SHADOW_TRACE_ABS_ACC = 1e-12;

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
            main_metric.calcThetaPhi(8.18 * unit.pc, theta, sin_theta, cos_theta, alpha, beta, u,
                                     result.mutable_data() + offset_0, result.mutable_data() + offset_1,
                                     result.mutable_data() + offset_2, result.mutable_data() + offset_3);
            offset_0 += offset;
            offset_1 += offset;
            offset_2 += offset;
            offset_3 += offset;
            bars[0].set_progress((i + 100 + (j + 100) / 201.) * 100. / 201.);
        }
    }
    return result.reshape({201, 201, 4, 201});
}

py::array_t<double> kerrRingShape(const py::array_t<double> r, double R, double spin, double theta) {
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
        main_metric.calcThetaPhi(8.18 * unit.pc, theta, sin_theta, cos_theta, alpha, beta, u, result.mutable_data() + offset_0,
                                 result.mutable_data() + offset_1, result.mutable_data() + offset_2,
                                 result.mutable_data() + offset_3);
        offset_0 += offset;
        offset_1 += offset;
        offset_2 += offset;
        offset_3 += offset;
    }
    return result.reshape({101, 4, N});
}

py::array_t<double> calculateShadowSilhouette(const py::array_t<double> black_hole_parameters, int metric) {
    return py::array_t<double>();
}

py::array_t<double> ShadowChi2(const py::array_t<double> observed_shadow_silhouette,
                               const py::array_t<double> black_hole_parameters, int metric, double resolution) {
    int data_point_number = observed_shadow_silhouette.shape(0);
    auto observed_shadow_silhouette_reader = observed_shadow_silhouette.unchecked<2>();
    vector<double> corresponding_weight(data_point_number);
    resolution = 0.5 / resolution;
    for (int i = 1; i < data_point_number; ++i) {
        double half_distance =
            resolution * sqrt(power2(observed_shadow_silhouette_reader(i, 0) - observed_shadow_silhouette_reader(i - 1, 0)) +
                              power2(observed_shadow_silhouette_reader(i, 1) - observed_shadow_silhouette_reader(i - 1, 1)));
        corresponding_weight[i] += half_distance;
        corresponding_weight[i - 1] += half_distance;
    }
    double half_distance =
        resolution *
        sqrt(power2(observed_shadow_silhouette_reader(data_point_number - 1, 0) - observed_shadow_silhouette_reader(0, 0)) +
             power2(observed_shadow_silhouette_reader(data_point_number - 1, 1) - observed_shadow_silhouette_reader(0, 1)));
    corresponding_weight[data_point_number - 1] += half_distance;
    corresponding_weight[0] += half_distance;
    auto result = py::array_t<double>(data_point_number);
    auto result_writer = result.mutable_unchecked();
    auto calculated_shadow_silhouette = calculateShadowSilhouette(black_hole_parameters, metric);
    auto calculated_shadow_silhouette_reader = calculated_shadow_silhouette.unchecked<2>();
    for (int i = 0; i < data_point_number; ++i) {
        result_writer(i) = corresponding_weight[i] *
                           sqrt(power2(calculated_shadow_silhouette_reader(i, 0) - observed_shadow_silhouette_reader(i, 0)) +
                                power2(calculated_shadow_silhouette_reader(i, 1) - observed_shadow_silhouette_reader(i, 1)));
    }
    return result;
}

void bindShadowModule(py::module_ &m) {
    m.def("kerr_coverage", &KerrCoverage, R"pbdoc(
		KerrCoverage
	)pbdoc");

    m.def("kerr_ring_shape", &kerrRingShape, R"pbdoc(
		KerrRingShape
	)pbdoc");
}
