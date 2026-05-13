/**
 * @file config.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-15
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "metric/metrics.hpp"

namespace sbody {
template <typename T>
struct RuntimeConfig {
    T mass_ = 4.15e6;
    std::string metric_name_ = "Schwarzschild";
    MetricType metric_;
    MetricParameters<T> metric_params_;
    T inclination_ = 0.0;
    T position_angle_ = 0.0;
    double t_final_ = 3600.;
    int t_step_number_ = 10000;
    int t_calculation_ = 3600000;
    int ray_ = 1;
    T absolute_accuracy_ = 1e-15;
    T relative_accuracy_ = 1e-15;
    std::string output_dir_ = "data";
    std::string output_format_ = "NumPy";
    std::vector<T> kepler_ = {1, 2, 3, 4, 5, 6};
    bool progress_bar_ = false;
    bool validate() {
        try {
            metric_ = getMetricType(metric_name_);
        } catch (const std::invalid_argument &e) {
            return false;
        }
        return true;
    }
};
} // namespace sbody
